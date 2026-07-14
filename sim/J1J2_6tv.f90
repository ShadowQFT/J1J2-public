program J1J2_6
   USE, INTRINSIC :: IEEE_EXCEPTIONS
   USE, INTRINSIC :: IEEE_ARITHMETIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Monte Carlo simulation of diluted 3d J1J2 model
!     using Metropolis updates
!
!     periodic boundary conditions
!     uses KISS05 random number generator
!     -------------------------------------------------------------------
!     History:
!
!     J1J2_1     :    6 Mar 2018  first version based on Ising2d9
!     J1J2_2     :    6 Mar 2018  first MPI version, conditional compilation
!     J1J2_3     :    8 Mar 2018  random and anticorrelated impurities
!     J1J2_4     :    1 Oct 2019  add stripe measurement, anisotropic J1 by Ye
!     J1J2_5     :   14 Mar 2021  add cornersweeps, nematic OP,
!                                 output of local config, domain energy
!     J1J2_6     :   21 Aug 2021  error bars for stripe OP + stripe susc.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
   integer, parameter    :: r8b = SELECTED_REAL_KIND(P=14, R=99)   ! 8-byte reals
   integer, parameter    :: i4b = SELECTED_INT_KIND(8)            ! 4-byte integers
   integer, parameter    :: i1b = SELECTED_INT_KIND(2)            ! 1-byte integers
   integer, parameter    :: ilog = kind(.true.)

   type observable_state
      real(r8b) :: magnetization
      real(r8b) :: energy
      real(r8b) :: stripe_x
      real(r8b) :: stripe_y
   end type observable_state

   type config_observable_summary
      real(r8b) :: mean_abs_magnetization
      real(r8b) :: magnetization_second_moment
      real(r8b) :: magnetization_fourth_moment
      real(r8b) :: mean_energy
      real(r8b) :: energy_second_moment
      real(r8b) :: mean_abs_stripe_x
      real(r8b) :: mean_abs_stripe_y
      real(r8b) :: stripe_x_second_moment
      real(r8b) :: stripe_x_fourth_moment
      real(r8b) :: stripe_y_second_moment
      real(r8b) :: stripe_y_fourth_moment
      real(r8b) :: stripe_second_moment
      real(r8b) :: stripe_fourth_moment
      real(r8b) :: mean_abs_stripe_order
      real(r8b) :: mean_nematic_order
      real(r8b) :: mean_local_nematic
      real(r8b) :: mean_abs_local_nematic
      real(r8b) :: local_nematic_second_moment
      real(r8b) :: local_nematic_fourth_moment
      real(r8b) :: best_energy
      real(r8b) :: final_energy
   end type config_observable_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer(i4b), parameter    :: L = 10                               ! in-plane linear system size
   integer(i4b), parameter    :: H = 10                               ! number of z layers
   real(r8b), parameter    :: J1 = 1.0D0                           ! nn interaction
   real(r8b), parameter    :: deltaJ1 = 0.0D0                      ! difference of nn interaction
   real(r8b), parameter    :: J1h = J1 + deltaJ1                   ! horizontal J of nn interaction
   real(r8b), parameter    :: J1v = J1 - deltaJ1                   ! vertical J of nn interaction
   real(r8b), parameter    :: J2 = -1.0D0                  ! nnn interaction in terms of nn interaction
   real(r8b), parameter    :: J3 = 1.0D0                           ! interaction between layers
   real(r8b), parameter    :: phi = 2.0D0                         ! domain strength
   real(r8b), parameter    :: corr_length = 3.0D0                  ! correlation length of random field
   real(r8b), parameter    :: TMAX = 14.0D0                          ! temperature loop control - maximum
   real(r8b), parameter    :: TMIN = 8.0D0                          ! temperature loop control - minimum
   real(r8b), parameter    :: DT = 2.0D0                          ! temperature step for LOCAL mode; PT derives spacing from TMIN/TMAX and NREPLICA
   real(r8b), parameter    :: PP = 0.0D0                          ! impurity concentration
   real(r8b), parameter    :: PPh = 0.0D0                            ! horizontal J impurity concentration
   real(r8b), parameter    :: PPv = 0.0D0                            ! vertical J impurity concentration
   character, parameter    :: IMPMODE*3 = 'RAN'                      ! RAN = random, COR = anticorrelated
   character, parameter    :: START*3 = 'STR'                        ! HOT or FER(ro) or STR(ipe)
   character, parameter    :: STRAXIS*1 = 'X'                        ! stripe axis: X or Y (stripe lines run parallel to this axis)
   character, parameter    :: UPDATEORDER*3 = 'SEQ'                  ! SEQ = fixed site order, RND = fresh random permutation each sweep
   character, parameter    :: TEMPPROTOCOL*5 = 'CARRY'               ! CARRY = reuse previous-T state, INDEP = fresh start at each temperature
   character, parameter    :: SAMPLERMODE*5 = 'LOCAL'                ! LOCAL = single-replica baseline, PT = parallel tempering
   character, parameter    :: PTSTARTMODE*6 = 'COMMON'               ! COMMON = shared start family, INDEP = per-replica reseeded starts
   character, parameter    :: PTLADDERMODE*6 = 'LINEAR'              ! LINEAR = uniform spacing, SPLIT = piecewise spacing around TLADDERMID
   character, parameter    :: PTPARALLELMODE*7 = 'SERIAL'            ! SERIAL = one rank owns full PT ladder, GROUPED = one rank per replica
   real(r8b), parameter    :: TLADDERMID = 9.5D0                     ! split-point temperature for nonlinear ladder
   real(r8b), parameter    :: TLADDERFRAC = 0.7D0                    ! fraction of ladder intervals assigned to [TMIN, TLADDERMID]
   integer(i4b), parameter    :: NEQ = 1000000                            ! Monte Carlo equilibration sweeps
   integer(i4b), parameter    :: NMESS = 1000000                          ! Monte Carlo measurement sweeps
   integer(i4b), parameter    :: NCONF = 100                           ! number of disorder configs
   integer(i4b), parameter    :: PTSWAPEVER = 5                       ! local sweeps between adjacent-replica swap passes in PT mode
   logical(ilog), parameter    :: WRITELOCNEM = .true.                  ! write file of local nematic OPs
   logical(ilog), parameter    :: WRITEHIST = .true.                    ! write histogram of energy
   logical(ilog), parameter    :: WRITESTRIP = .true.                ! write strip x and strip y for each sample
   logical(ilog), parameter    :: WRITENUMFLIPS = .false.                ! write numflips file
   logical(ilog), parameter    :: WRITEBINDTIME = .true.                 ! write low-T binder traces for selected disorder configs
   logical(ilog), parameter    :: WRITEPHITABLE = .true.                 ! write phi at each site
   logical(ilog), parameter    :: MAPTOUNIFORM = .true.                  ! map random field to uniform distribution
   logical(ilog), parameter    :: VALIDATE_INCREMENTAL = .false.         ! cross-check maintained observables against full measure_observables

   integer(i4b), parameter    :: BINS = 100                             ! number of bins for histogram
   integer(i4b), parameter    :: NTRACECONF = 5                         ! number of disorder configs for low-T binder trace output

   integer(i4b), parameter    :: IRINIT = 1                             ! Random number seed, must be positive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants - do not touch !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer(i4b), parameter    :: L2 = L*L                                    ! number of sites in one layer
   integer(i4b), parameter    :: L3 = L2*H                                   ! total number of sites
   integer(i4b), parameter    :: NREPLICA = 1                                ! PT uses one replica per ladder temperature
   integer(i4b), parameter    :: NTEMP = merge(NREPLICA, 1 + NINT((TMAX - TMIN)/DT), SAMPLERMODE == 'PT')
   integer(i4b), parameter    :: PTMEASINDEX = 1                             ! target temperature index for PT diagnostics (1 = TMIN)

   type temperature_measurement_sums
      real(r8b) :: abs_magnetization(NTEMP)
      real(r8b) :: magnetization_second_moment(NTEMP)
      real(r8b) :: magnetization_fourth_moment(NTEMP)
      real(r8b) :: energy(NTEMP)
      real(r8b) :: energy_second_moment(NTEMP)
      real(r8b) :: abs_stripe_x(NTEMP)
      real(r8b) :: abs_stripe_y(NTEMP)
      real(r8b) :: stripe_x_second_moment(NTEMP)
      real(r8b) :: stripe_x_fourth_moment(NTEMP)
      real(r8b) :: stripe_y_second_moment(NTEMP)
      real(r8b) :: stripe_y_fourth_moment(NTEMP)
      real(r8b) :: stripe_second_moment(NTEMP)
      real(r8b) :: stripe_fourth_moment(NTEMP)
      real(r8b) :: abs_stripe_order(NTEMP)
      real(r8b) :: nematic_order(NTEMP)
      real(r8b) :: local_nematic(NTEMP)
      real(r8b) :: abs_local_nematic(NTEMP)
      real(r8b) :: local_nematic_second_moment(NTEMP)
      real(r8b) :: local_nematic_fourth_moment(NTEMP)
      real(r8b) :: best_energy(NTEMP)
      real(r8b) :: final_energy(NTEMP)
   end type temperature_measurement_sums

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variable declarations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Naming rule used below:
! sample_*   = one measurement sweep inside one temperature/config
! config_*   = averages or diagnostics for the current disorder config
! disorder_* = sums over disorder configs for the final output files

   real(r8b)      :: stripe_by_config(NCONF, 2, NTEMP)                        ! stripex and stripey for each disorder config
   real(r8b)      :: config_stripe_x(NTEMP), config_stripe_y(NTEMP)                  ! current-config stripe averages by temperature
   real(r8b)      :: config_abs_stripe(NTEMP), config_stripe_difference(NTEMP)     ! per-config stripe diagnostics by temperature
   real(r8b)      :: config_best_energy(NTEMP), config_final_energy(NTEMP)             ! per-config best/final energies by temperature
   real(r8b)      :: abs_stripe_by_config(NCONF, NTEMP), stripe_difference_by_config(NCONF, NTEMP)
   real(r8b)      :: best_energy_by_config(NCONF, NTEMP), final_energy_by_config(NCONF, NTEMP)

   integer(i4b)    :: flip_counts_this_config(NTEMP, 0:L3 - 1)   ! number of flips per site at each temperature
   integer(i4b)    :: flip_counts_by_config(NCONF, NTEMP, 0:L3 - 1)    ! number of flips per site in each disorder config
   integer(i4b)    :: config_flip_counts(NTEMP, 0:L3 - 1)             ! current-config flip counts by temperature

   real(r8b)       :: config_energy_histogram(NTEMP, BINS)                           ! list of energy histogram values
   real(r8b)       :: disorder_energy_histogram(NTEMP, BINS)
   real(r8b)       :: energy_histogram_min(NTEMP), energy_histogram_max(NTEMP), energy_histogram_bin_size(NTEMP)

   integer(i1b)    :: s1(0:L3 - 1)                                      ! Ising variable
   integer(i1b)    :: s1_rep(0:L3 - 1, NREPLICA)                         ! PT replica spin states
   integer(i1b)    :: jh(0:L3 - 1)                                      ! horizontal J validity
   integer(i1b)    :: jv(0:L3 - 1)                                      ! vertical J validity
   integer(i1b)    :: occu(0:L3 - 1)                                    ! occu=1 site occupied, occu=0 site empty
   real(r8b)       :: sample_magnetization(NMESS), sample_energy(NMESS)                    ! list of magnetization and energy for each sweep
   real(r8b)       :: sample_stripe_x(NMESS), sample_stripe_y(NMESS)                 ! list of stripex and stripey for each sweep
   real(r8b)       :: sample_local_nematic_mean(NMESS)                                  ! list of spatially averaged local nematic per sweep

   real(r8b)       :: config_abs_magnetization(NTEMP)       ! current-config magnetization averages
   real(r8b)       :: config_abs_magnetization_squared(NTEMP), config_log_abs_magnetization(NTEMP)
   real(r8b)       :: config_magnetization_second_moment(NTEMP)
   real(r8b)       :: config_magnetization_fourth_moment(NTEMP)
   real(r8b)       :: config_abs_stripe_order(NTEMP), config_abs_stripe_order_squared(NTEMP)
   real(r8b)       :: config_stripe_second_moment(NTEMP), config_stripe_fourth_moment(NTEMP)
   real(r8b)       :: config_abs_stripe_x(NTEMP), config_abs_stripe_y(NTEMP)
   real(r8b)       :: config_stripe_x_second_moment(NTEMP), config_stripe_y_second_moment(NTEMP)
   real(r8b)       :: config_stripe_x_fourth_moment(NTEMP), config_stripe_y_fourth_moment(NTEMP)
   real(r8b)       :: config_magnetic_susceptibility(NTEMP), config_magnetic_susceptibility_squared(NTEMP)
   real(r8b)       :: config_stripe_susceptibility(NTEMP), config_stripe_susceptibility_squared(NTEMP)
   real(r8b)       :: config_magnetization_binder(NTEMP), config_magnetization_binder_squared(NTEMP)
   real(r8b)       :: config_stripe_binder(NTEMP), config_stripe_binder_squared(NTEMP)
   real(r8b)       :: config_stripe_x_binder(NTEMP), config_stripe_x_binder_squared(NTEMP)
   real(r8b)       :: config_stripe_y_binder(NTEMP), config_stripe_y_binder_squared(NTEMP)
   real(r8b)       :: config_energy(NTEMP), config_energy_squared(NTEMP)                   ! configuration averages of energy
   real(r8b)       :: config_specific_heat(NTEMP), config_specific_heat_squared(NTEMP)                  ! configuration averages of specific heat
   real(r8b)       :: config_nematic_order(NTEMP), config_nematic_order_squared(NTEMP)                 ! configuration averages of nematic OP
   real(r8b)       :: config_local_nematic_mean(NTEMP), config_abs_local_nematic_mean(NTEMP)
   real(r8b)       :: config_local_nematic_mean_squared(NTEMP)
   real(r8b)       :: config_local_nematic_second_moment(NTEMP), config_local_nematic_fourth_moment(NTEMP)
   real(r8b)       :: config_local_nematic_binder(NTEMP), config_local_nematic_binder_squared(NTEMP)

   real(r8b)       :: disorder_abs_magnetization(NTEMP)       ! sums over disorder configs
   real(r8b)       :: disorder_abs_magnetization_squared(NTEMP), disorder_log_abs_magnetization(NTEMP)
   real(r8b)       :: disorder_magnetization_second_moment(NTEMP)
   real(r8b)       :: disorder_magnetization_fourth_moment(NTEMP)
   real(r8b)       :: disorder_abs_stripe_order(NTEMP), disorder_abs_stripe_order_squared(NTEMP)
   real(r8b)       :: disorder_stripe_second_moment(NTEMP), disorder_stripe_fourth_moment(NTEMP)
   real(r8b)       :: disorder_abs_stripe_x(NTEMP), disorder_abs_stripe_y(NTEMP)
   real(r8b)       :: disorder_stripe_x_second_moment(NTEMP), disorder_stripe_y_second_moment(NTEMP)
   real(r8b)       :: disorder_stripe_x_fourth_moment(NTEMP), disorder_stripe_y_fourth_moment(NTEMP)
   real(r8b)       :: disorder_magnetic_susceptibility(NTEMP), disorder_magnetic_susceptibility_squared(NTEMP)
   real(r8b)       :: disorder_stripe_susceptibility(NTEMP), disorder_stripe_susceptibility_squared(NTEMP)
   real(r8b)       :: disorder_magnetization_binder(NTEMP), disorder_magnetization_binder_squared(NTEMP)
   real(r8b)       :: disorder_stripe_binder(NTEMP), disorder_stripe_binder_squared(NTEMP)
   real(r8b)       :: disorder_stripe_x_binder(NTEMP), disorder_stripe_x_binder_squared(NTEMP)
   real(r8b)       :: disorder_stripe_y_binder(NTEMP), disorder_stripe_y_binder_squared(NTEMP)
   real(r8b)       :: disorder_energy(NTEMP), disorder_energy_squared(NTEMP)                     ! configuration averages of energy
   real(r8b)       :: disorder_specific_heat(NTEMP), disorder_specific_heat_squared(NTEMP)                    ! configuration averages of specific heat
   real(r8b)       :: disorder_nematic_order(NTEMP), disorder_nematic_order_squared(NTEMP)
   real(r8b)       :: disorder_local_nematic_mean(NTEMP), disorder_abs_local_nematic_mean(NTEMP)
   real(r8b)       :: disorder_local_nematic_mean_squared(NTEMP)
   real(r8b)       :: disorder_local_nematic_second_moment(NTEMP), disorder_local_nematic_fourth_moment(NTEMP)
   real(r8b)       :: disorder_local_nematic_binder(NTEMP), disorder_local_nematic_binder_squared(NTEMP)

   real(r8b)       :: local_nematic_mean_by_config(NCONF, NTEMP)
   real(r8b)       :: abs_local_nematic_mean_by_config(NCONF, NTEMP), local_nematic_second_moment_by_config(NCONF, NTEMP)
   real(r8b)       :: local_nematic_fourth_moment_by_config(NCONF, NTEMP), local_nematic_binder_by_config(NCONF, NTEMP)

   real(r8b)       :: p(0:L3 - 1)             ! random field at each site
   real(r8b)       :: phi_by_config(NCONF, 0:L3 - 1)        ! joined list of random fields over all disorder configs
   real(r8b)       :: config_phi(0:L3 - 1)               ! current-config random field values

   integer(i4b)    :: m1(0:L3 - 1)            ! neighbor table: nearest neighbors
   integer(i4b)    :: m2(0:L3 - 1)
   integer(i4b)    :: m3(0:L3 - 1)
   integer(i4b)    :: m4(0:L3 - 1)
   integer(i4b)    :: m5(0:L3 - 1)            ! neighbor table: next nearest neighbors
   integer(i4b)    :: m6(0:L3 - 1)
   integer(i4b)    :: m7(0:L3 - 1)
   integer(i4b)    :: m8(0:L3 - 1)
   integer(i4b)    :: m9(0:L3 - 1)            ! neighbor table: next layer neighbors
   integer(i4b)    :: m10(0:L3 - 1)

   integer(i4b)    :: iconf                 ! current disorder config
   integer(i4b)    :: confstep             ! current disorder-loop step under MPI
   integer(i4b)    :: max_cfg_per_rank     ! maximum number of disorder configs assigned to any rank
   integer(i4b)    :: completed_config_count
   integer(i4b)    :: ibin                  ! current bin
   integer(i4b)    :: nimp                  ! number of impurity sites
   integer(i4b)    :: isite                 ! current site index
   integer(i4b)    :: layer                 !  layer index

   real(r8b)       :: T, beta               ! Monte Carlo temperature
   real(r8b)       :: temp_ladder(NTEMP), beta_ladder(NTEMP)            ! explicit ladder used by PT and output helpers
   real(r8b)       :: spin_sum_by_site(0:L3 - 1), local_nematic_sum_by_site(0:L3 - 1)   ! measurement sums for site diagnostics
   real(r8b)       :: stripe_parity_x(0:L3 - 1), stripe_parity_y(0:L3 - 1)
   real(r8b)       :: config_local_nematic_by_temp(NTEMP), config_abs_local_nematic_by_temp(NTEMP)
   real(r8b)       :: config_local_nematic_second_by_temp(NTEMP), config_local_nematic_fourth_by_temp(NTEMP)
   real(r8b)       :: config_local_nematic_binder_by_temp(NTEMP)
   integer(i4b)    :: local_nematic_sample_count
   integer(i4b)    :: temperature_index
   integer(i4b)    :: site_order(0:L3 - 1)
   integer(i4b)    :: replica_at_temperature(NREPLICA), temperature_of_replica(NREPLICA)      ! PT ownership maps between configs and ladder slots
   integer(i4b)    :: config_pt_swap_attempts(NREPLICA - 1), config_pt_swap_accepts(NREPLICA - 1)
   integer(i4b)    :: disorder_pt_swap_attempts(NREPLICA - 1), disorder_pt_swap_accepts(NREPLICA - 1)
   integer(i4b)    :: config_pt_visits_both_ends, disorder_pt_visits_both_ends

   integer(i4b)    :: init
   integer(kind=8) :: config_clock_start, collect_clock_start
   logical(ilog)   :: active_config

   real(r8b), external     :: rkiss05
   external kissinit

! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
   include 'mpif.h'
   integer(i4b)              :: ierr
   integer(i4b)              :: source_rank, myid         ! CPU index
   integer(i4b)              :: numprocs                 ! total number of CPUs
   integer(i4b)              :: status(MPI_STATUS_SIZE)
   integer(i4b)              :: pt_comm
   integer(i4b)              :: pt_group_id, pt_group_rank, pt_group_size
   integer(i4b)              :: pt_num_groups, pt_active_ranks
   logical(ilog)             :: pt_group_active, pt_group_master
#endif

   character cfenfile*16, cfmagfile*16, cfstrfile*16, cfnmfile*16, cfhsfile*16, cfxyfile*16, cfndfile*16, cfnumfile*16, cfptable*16, cfbtfile*16, cfptfile*16
   cfenfile = 'cfen00000000.dat'
   cfmagfile = 'cfma00000000.dat'
   cfstrfile = 'cfst00000000.dat'
   cfnmfile = 'cfnm00000000.dat'
   cfhsfile = 'cfhs00000000.dat'
   cfxyfile = 'cfxy00000000.dat'
   cfndfile = 'cfnd00000000.dat'
   cfnumfile = 'cfnu00000000.dat'
   cfbtfile = 'cfbt00000000.dat'
   cfptfile = 'cfpt00000000.dat'
   cfptable = 'cfph00000000.dat'
   write (cfenfile(5:8), '(I4.4)') L
   write (cfenfile(9:12), '(I4.4)') H
   write (cfmagfile(5:8), '(I4.4)') L
   write (cfmagfile(9:12), '(I4.4)') H
   write (cfstrfile(5:8), '(I4.4)') L
   write (cfstrfile(9:12), '(I4.4)') H
   write (cfnmfile(5:8), '(I4.4)') L
   write (cfnmfile(9:12), '(I4.4)') H
   write (cfhsfile(5:8), '(I4.4)') L
   write (cfhsfile(9:12), '(I4.4)') H
   write (cfxyfile(5:8), '(I4.4)') L
   write (cfxyfile(9:12), '(I4.4)') H
   write (cfndfile(5:8), '(I4.4)') L
   write (cfndfile(9:12), '(I4.4)') H
   write (cfnumfile(5:8), '(I4.4)') L
   write (cfnumfile(9:12), '(I4.4)') H
   write (cfbtfile(5:8), '(I4.4)') L
   write (cfbtfile(9:12), '(I4.4)') H
   write (cfptfile(5:8), '(I4.4)') L
   write (cfptfile(9:12), '(I4.4)') H
   write (cfptable(5:8), '(I4.4)') L
   write (cfptable(9:12), '(I4.4)') H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start of main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call IEEE_SET_UNDERFLOW_MODE(.TRUE.)

#ifdef PARALLEL
   call setup_mpi                                  ! Start MPI
   call setup_pt_parallel
#endif

   call setup_neighbor_table                       ! build neighbor table for the 3D periodic lattice
   call setup_stripe_parities
   call setup_config_averages                        ! zero the arrays containing configuration averages
   call validate_run_options
   call setup_temperature_ladder

#ifdef PARALLEL
   call set_disorder_loop_length
   disorder_loop: do confstep = 0, max_cfg_per_rank - 1
      call assign_disorder_config                  ! pick which disorder config this rank/group owns
#else
      disorder_loop: do iconf = 1, NCONF
         call assign_disorder_config
#endif
         if (.not. active_config) then
#ifdef PARALLEL
            ! inactive ranks still join the collect call so MPI does not hang
            call collect_config_data(active_config)
#endif
            cycle disorder_loop
         end if
         call setup_disorder_config                 ! reset arrays, seed RNG, make vacancies and phi field
         call run_sampler_for_config                ! do either local Metropolis or PT for this config
         call collect_config_results                ! send/merge this config into the running sums
         call write_current_results                 ! rewrite output after each finished batch

      end do disorder_loop

#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      stop   ! Main Program stops !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Now the internal subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      subroutine set_disorder_loop_length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! grouped PT uses one MPI group per disorder config, not one rank
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
            max_cfg_per_rank = (NCONF + pt_num_groups - 1)/pt_num_groups
         else
            max_cfg_per_rank = (NCONF + numprocs - 1)/numprocs
         end if
      end subroutine set_disorder_loop_length
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine assign_disorder_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
         ! In parallel, iconf is assigned here instead of by the do-loop index.
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
            if (pt_group_active) then
               iconf = pt_group_id + 1 + confstep*pt_num_groups
            else
               iconf = 0
            end if
            active_config = (pt_group_active .and. iconf <= NCONF)
            if (pt_group_master .and. active_config) print *, 'dis. conf.', iconf
         else
            iconf = myid + 1 + confstep*numprocs
            active_config = (iconf <= NCONF)
            if (myid == 0 .and. active_config) print *, 'dis. conf.', iconf
         end if
#else
         ! Serial mode just runs configs in order.
         active_config = .true.
         print *, 'dis. conf.', iconf
#endif
      end subroutine assign_disorder_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_disorder_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call reset_current_config_averages           ! zero the arrays containing current-config averages
         call timer_start(config_clock_start)         ! start timing for this disorder config
         call write_progress('CONFIG_START', 0, 0)
         init = IRINIT + iconf - 1
         call kissinit(init)                          ! same config number gives same random start
         call setup_vacancies                         ! draw missing sites for this disorder config
         call mark_vacancy_flipcounts
         if (SAMPLERMODE == 'LOCAL' .and. TEMPPROTOCOL == 'CARRY') then
            call initialize_spin_array(s1)             ! CARRY starts once and reuses state across T
         end if
         call setup_phi_table                         ! make the random field used by the Hamiltonian
      end subroutine setup_disorder_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine run_sampler_for_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Keep the mode decision here so the main loop stays readable.
         if (SAMPLERMODE == 'PT') then
            call run_parallel_tempering_config
         else
            call run_local_config
         end if
      end subroutine run_sampler_for_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine run_parallel_tempering_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! PT does the whole temperature ladder at once for this disorder config.
         if (PTPARALLELMODE == 'GROUPED') then
#ifdef PARALLEL
            call grouped_pt_run
            if (WRITELOCNEM .eqv. .true. .and. iconf == 1 .and. pt_group_master) call write_local_nematic_order
#else
            print *, 'ERROR: PTPARALLELMODE=GROUPED requires PARALLEL build.'
            stop 1
#endif
         else
            call initialize_pt_state                  ! one replica per ladder temperature
            call pt_equilibration_loop                ! local sweeps plus swap attempts
            call pt_measurement_loop                  ! same, but accumulate observables
            if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call snapshot_pt_target_replica
            if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call write_local_nematic_order
         end if
      end subroutine run_parallel_tempering_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine run_local_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Local mode still steps through the temperatures one at a time.
         do temperature_index = 1, NTEMP
            call run_local_temperature
         end do
      end subroutine run_local_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine run_local_temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         T = temperature_for_index(temperature_index)
         if (myid == 0) then
            print *, 'Running T= ', T, 'temperature_index= ', temperature_index
         end if
         call prepare_temperature_state
         beta = 1.D0/T
         flip_counts_this_config(temperature_index, :) = 0                  ! reset flip counters for this temperature
         where (occu(:) == 0)
            flip_counts_this_config(temperature_index, :) = -1              ! vacancies are marked so they are easy to skip later
         end where

         call equilibration_loop                      ! throw away sweeps before measurements
         call measurement_loop                        ! save per-sweep mag, energy, stripe, nematic
         call compute_config_observables              ! turn per-sweep lists into config averages

         if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call write_local_nematic_order
      end subroutine run_local_temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine collect_config_results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call timer_start(collect_clock_start)
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
#ifdef PARALLEL
            ! Only the group master owns a complete grouped-PT result to send up.
            call collect_config_data(active_config .and. pt_group_master)
#else
            call collect_config_data(active_config)
#endif
         else
            call collect_config_data(active_config)
         end if
         call write_timing('COLLECT_DATA', collect_clock_start)
         call write_timing('CONFIG_TOTAL', config_clock_start)
      end subroutine collect_config_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_current_results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
         ! These files are partial while the run is still going, but useful if a job dies.
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
            completed_config_count = min(NCONF, (confstep + 1)*pt_num_groups)
         else
            completed_config_count = min(NCONF, (confstep + 1)*numprocs)
         end if
#else
         completed_config_count = iconf
#endif
         call write_output_files(completed_config_count)
      end subroutine write_current_results

#ifdef PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

      end subroutine setup_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_pt_parallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: color

         pt_comm = MPI_COMM_NULL
         pt_group_id = myid
         pt_group_rank = 0
         pt_group_size = 1
         pt_num_groups = numprocs
         pt_active_ranks = numprocs
         pt_group_active = .true.
         pt_group_master = .true.

         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
            pt_num_groups = numprocs/NREPLICA
            pt_active_ranks = pt_num_groups*NREPLICA
            pt_group_active = (myid < pt_active_ranks)
            if (pt_num_groups < 1) then
               if (myid == 0) print *, 'ERROR: GROUPED PT requires at least NREPLICA MPI ranks.'
               call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            end if
            if (pt_group_active) then
               pt_group_id = myid/NREPLICA
               pt_group_rank = mod(myid, NREPLICA)
               color = pt_group_id
            else
               pt_group_id = -1
               pt_group_rank = -1
               color = MPI_UNDEFINED
            end if
            pt_group_master = (pt_group_active .and. pt_group_rank == 0)
            call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, myid, pt_comm, ierr)
            if (pt_group_active) then
               call MPI_COMM_RANK(pt_comm, pt_group_rank, ierr)
               call MPI_COMM_SIZE(pt_comm, pt_group_size, ierr)
            end if
         end if

      end subroutine setup_pt_parallel
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine timer_start(timer_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(kind=8), intent(out) :: timer_count

         call system_clock(timer_count)
      end subroutine timer_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_timing(label, timer_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         character(len=*), intent(in) :: label
         integer(kind=8), intent(in)  :: timer_count
         integer(kind=8)              :: now_count, count_rate
         real(r8b)                    :: elapsed

         call system_clock(now_count, count_rate)
         if (count_rate > 0) then
            elapsed = real(now_count - timer_count, r8b)/real(count_rate, r8b)
         else
            elapsed = -1.D0
         end if
#ifdef PARALLEL
         write (*, '(A,1x,A,1x,A,I0,1x,A,I0,1x,A,F12.3)') 'timing', trim(label), 'rank', myid, &
     &      'iconf', iconf, 'seconds', elapsed
#else
         write (*, '(A,1x,A,1x,A,I0,1x,A,I0,1x,A,F12.3)') 'timing', trim(label), 'rank', 0, &
     &      'iconf', iconf, 'seconds', elapsed
#endif
         call flush (6)
      end subroutine write_timing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_progress(phase, current, total)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         character(len=*), intent(in) :: phase
         integer(i4b), intent(in)     :: current, total
         character(len=64)            :: job_id

#ifdef PARALLEL
         if (myid /= 0) return
#endif

         job_id = ''
         call get_environment_variable('SLURM_JOB_ID', job_id)
         open (98, file='progress.txt', status='replace')
         write (98, '(A,1x,A)') 'job_id', trim(job_id)
         write (98, '(A,1x,A)') 'sampler', trim(SAMPLERMODE)
         write (98, '(A,1x,A)') 'phase', trim(phase)
         write (98, '(A,1x,I0)') 'iconf', iconf
         write (98, '(A,1x,I0)') 'nconf', NCONF
#ifdef PARALLEL
         write (98, '(A,1x,I0)') 'ntasks', numprocs
#else
         write (98, '(A,1x,I0)') 'ntasks', 1
#endif
         write (98, '(A,1x,I0)') 'current', current
         write (98, '(A,1x,I0)') 'total', total
         close (98)

         write (*, '(A,1x,A,1x,A,I0,1x,A,I0,1x,A,I0)') 'progress', trim(phase), 'iconf', iconf, &
     &      'step', current, 'of', total
         call flush (6)
      end subroutine write_progress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_neighbor_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)     :: ix, iy, iz                   ! x, y, and z coordinates

         do iz = 0, H - 1
            do iy = 0, L - 1
               do ix = 0, L - 1
                  isite = L*iy + ix + L2*iz
                  if (iy .eq. L - 1) then
                     m1(isite) = isite - L*(L - 1)
                  else
                     m1(isite) = isite + L
                  end if

                  if (iy .eq. 0) then
                     m2(isite) = isite + L*(L - 1)
                  else
                     m2(isite) = isite - L
                  end if

                  if (ix .eq. L - 1) then
                     m3(isite) = isite - (L - 1)
                  else
                     m3(isite) = isite + 1
                  end if

                  if (ix .eq. 0) then
                     m4(isite) = isite + (L - 1)
                  else
                     m4(isite) = isite - 1
                  end if

                  if (iz .eq. (H - 1)) then
                     m9(isite) = isite - L2*(H - 1)
                  else
                     m9(isite) = isite + L2
                  end if

                  if (iz .eq. 0) then
                     m10(isite) = isite + L2*(H - 1)
                  else
                     m10(isite) = isite - L2
                  end if
               end do
            end do
         end do

         do isite = 0, L3 - 1
            m5(isite) = m3(m1(isite))
            m6(isite) = m4(m1(isite))
            m7(isite) = m3(m2(isite))
            m8(isite) = m4(m2(isite))
         end do

      end subroutine setup_neighbor_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_phi_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         use corr_rand_3d_mod, only: make_corr_rand_3d
         implicit none

         integer, parameter :: r8b = selected_real_kind(14, 99)
         integer, parameter :: i4b = selected_int_kind(8)

         integer(i4b) :: ix, iy, layer

         integer :: n1, n2, n3
         real(r8b), allocatable :: field(:, :, :)
         real(r8b) :: xi
         integer(i4b) :: seed
         logical :: do_uniform
         real(r8b) :: max_value

         n1 = L; n2 = L; n3 = H
         allocate (field(0:n1 - 1, 0:n2 - 1, 0:n3 - 1))

         xi = corr_length
         seed = IRINIT + iconf - 1
         do_uniform = MAPTOUNIFORM
         max_value = phi

         call make_corr_rand_3d(field, n1, n2, n3, xi, seed, do_uniform, max_value)

         do layer = 0, H - 1
            do iy = 0, L - 1
               do ix = 0, L - 1
                  isite = L*iy + ix + L2*layer
                  p(isite) = field(ix, iy, layer)
               end do
            end do
         end do

         close (99)

      end subroutine setup_phi_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_config_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         stripe_by_config(:, :, :) = 0.D0
         config_stripe_x(:) = 0.D0
         config_stripe_y(:) = 0.D0
         config_abs_stripe(:) = 0.D0
         config_stripe_difference(:) = 0.D0
         config_best_energy(:) = 0.D0
         config_final_energy(:) = 0.D0
         config_local_nematic_by_temp(:) = 0.D0
         config_abs_local_nematic_by_temp(:) = 0.D0
         config_local_nematic_second_by_temp(:) = 0.D0
         config_local_nematic_fourth_by_temp(:) = 0.D0
         config_local_nematic_binder_by_temp(:) = 0.D0
         abs_stripe_by_config(:, :) = 0.D0
         stripe_difference_by_config(:, :) = 0.D0
         best_energy_by_config(:, :) = 0.D0
         final_energy_by_config(:, :) = 0.D0
         local_nematic_mean_by_config(:, :) = 0.D0
         abs_local_nematic_mean_by_config(:, :) = 0.D0
         local_nematic_second_moment_by_config(:, :) = 0.D0
         local_nematic_fourth_moment_by_config(:, :) = 0.D0
         local_nematic_binder_by_config(:, :) = 0.D0

         flip_counts_this_config(:, :) = 0
         flip_counts_by_config(:, :, :) = 0
         config_flip_counts(:, :) = 0

         config_abs_magnetization(:) = 0.D0
         config_magnetization_second_moment(:) = 0.D0
         config_magnetization_fourth_moment(:) = 0.D0
         config_magnetic_susceptibility(:) = 0.D0
         config_magnetization_binder(:) = 0.D0
         config_magnetic_susceptibility_squared(:) = 0.D0
         config_magnetization_binder_squared(:) = 0.D0
         config_energy(:) = 0.D0
         config_specific_heat(:) = 0.D0
         config_energy_squared(:) = 0.D0
         config_specific_heat_squared(:) = 0.D0
         config_stripe_binder(:) = 0.D0
         config_stripe_binder_squared(:) = 0.D0
         config_stripe_second_moment(:) = 0.D0
         config_stripe_fourth_moment(:) = 0.D0
         config_abs_stripe_order(:) = 0.D0
         config_abs_stripe_x(:) = 0.D0
         config_abs_stripe_y(:) = 0.D0
         config_stripe_x_second_moment(:) = 0.D0
         config_stripe_y_second_moment(:) = 0.D0
         config_stripe_x_fourth_moment(:) = 0.D0
         config_stripe_y_fourth_moment(:) = 0.D0
         config_abs_stripe_order_squared(:) = 0.D0
         config_stripe_susceptibility(:) = 0.D0
         config_stripe_susceptibility_squared(:) = 0.D0
         config_log_abs_magnetization(:) = 0.D0
         config_nematic_order(:) = 0.D0
         config_nematic_order_squared(:) = 0.D0
         config_stripe_x_binder(:) = 0.D0
         config_stripe_x_binder_squared(:) = 0.D0
         config_stripe_y_binder(:) = 0.D0
         config_stripe_y_binder_squared(:) = 0.D0
         config_local_nematic_mean(:) = 0.D0
         config_abs_local_nematic_mean(:) = 0.D0
         config_local_nematic_mean_squared(:) = 0.D0
         config_local_nematic_second_moment(:) = 0.D0
         config_local_nematic_fourth_moment(:) = 0.D0
         config_local_nematic_binder(:) = 0.D0
         config_local_nematic_binder_squared(:) = 0.D0

         disorder_abs_magnetization(:) = 0.D0
         disorder_magnetization_second_moment(:) = 0.D0
         disorder_magnetization_fourth_moment(:) = 0.D0
         disorder_magnetic_susceptibility(:) = 0.D0
         disorder_magnetic_susceptibility_squared(:) = 0.D0
         disorder_magnetization_binder(:) = 0.D0
         disorder_magnetization_binder_squared(:) = 0.D0
         disorder_abs_magnetization_squared(:) = 0.D0
         disorder_log_abs_magnetization(:) = 0.D0
         disorder_energy(:) = 0.D0
         disorder_energy_squared(:) = 0.D0
         disorder_specific_heat(:) = 0.D0
         disorder_specific_heat_squared(:) = 0.D0
         disorder_abs_stripe_order(:) = 0.D0
         disorder_abs_stripe_x(:) = 0.D0
         disorder_abs_stripe_y(:) = 0.D0
         disorder_stripe_second_moment(:) = 0.D0
         disorder_stripe_fourth_moment(:) = 0.D0
         disorder_stripe_susceptibility(:) = 0.D0
         disorder_stripe_susceptibility_squared(:) = 0.D0
         disorder_stripe_binder(:) = 0.D0
         disorder_stripe_binder_squared(:) = 0.D0
         disorder_nematic_order(:) = 0.D0
         disorder_nematic_order_squared(:) = 0.D0
         disorder_abs_stripe_order_squared(:) = 0.D0
         disorder_stripe_x_second_moment(:) = 0.D0
         disorder_stripe_x_fourth_moment(:) = 0.D0
         disorder_stripe_x_binder(:) = 0.D0
         disorder_stripe_x_binder_squared(:) = 0.D0
         disorder_stripe_y_second_moment(:) = 0.D0
         disorder_stripe_y_fourth_moment(:) = 0.D0
         disorder_stripe_y_binder(:) = 0.D0
         disorder_stripe_y_binder_squared(:) = 0.D0
         disorder_local_nematic_mean(:) = 0.D0
         disorder_abs_local_nematic_mean(:) = 0.D0
         disorder_local_nematic_mean_squared(:) = 0.D0
         disorder_local_nematic_second_moment(:) = 0.D0
         disorder_local_nematic_fourth_moment(:) = 0.D0
         disorder_local_nematic_binder(:) = 0.D0
         disorder_local_nematic_binder_squared(:) = 0.D0

         config_energy_histogram(:, :) = 0.D0
         disorder_energy_histogram(:, :) = 0.D0
         energy_histogram_min(:) = 0.D0
         energy_histogram_max(:) = 0.D0
         energy_histogram_bin_size(:) = 0.D0
         config_phi(:) = 0.D0
         s1_rep(:, :) = 1
         replica_at_temperature(:) = 0
         temperature_of_replica(:) = 0
         config_pt_swap_attempts(:) = 0
         config_pt_swap_accepts(:) = 0
         disorder_pt_swap_attempts(:) = 0
         disorder_pt_swap_accepts(:) = 0
         config_pt_visits_both_ends = 0
         disorder_pt_visits_both_ends = 0
         temp_ladder(:) = 0.D0
         beta_ladder(:) = 0.D0
         spin_sum_by_site(:) = 0.D0
         local_nematic_sum_by_site(:) = 0.D0
         local_nematic_sample_count = 0

      end subroutine setup_config_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reset_current_config_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         config_stripe_x(:) = 0.D0
         config_stripe_y(:) = 0.D0
         config_abs_stripe(:) = 0.D0
         config_stripe_difference(:) = 0.D0
         config_best_energy(:) = 0.D0
         config_final_energy(:) = 0.D0
         config_local_nematic_by_temp(:) = 0.D0
         config_abs_local_nematic_by_temp(:) = 0.D0
         config_local_nematic_second_by_temp(:) = 0.D0
         config_local_nematic_fourth_by_temp(:) = 0.D0
         config_local_nematic_binder_by_temp(:) = 0.D0

         flip_counts_this_config(:, :) = 0
         config_flip_counts(:, :) = 0

         config_abs_magnetization(:) = 0.D0
         config_magnetization_second_moment(:) = 0.D0
         config_magnetization_fourth_moment(:) = 0.D0
         config_magnetic_susceptibility(:) = 0.D0
         config_magnetization_binder(:) = 0.D0
         config_magnetic_susceptibility_squared(:) = 0.D0
         config_magnetization_binder_squared(:) = 0.D0
         config_energy(:) = 0.D0
         config_specific_heat(:) = 0.D0
         config_energy_squared(:) = 0.D0
         config_specific_heat_squared(:) = 0.D0
         config_stripe_binder(:) = 0.D0
         config_stripe_binder_squared(:) = 0.D0
         config_stripe_second_moment(:) = 0.D0
         config_stripe_fourth_moment(:) = 0.D0
         config_abs_stripe_order(:) = 0.D0
         config_abs_stripe_x(:) = 0.D0
         config_abs_stripe_y(:) = 0.D0
         config_stripe_x_second_moment(:) = 0.D0
         config_stripe_y_second_moment(:) = 0.D0
         config_stripe_x_fourth_moment(:) = 0.D0
         config_stripe_y_fourth_moment(:) = 0.D0
         config_abs_stripe_order_squared(:) = 0.D0
         config_stripe_susceptibility(:) = 0.D0
         config_stripe_susceptibility_squared(:) = 0.D0
         config_log_abs_magnetization(:) = 0.D0
         config_nematic_order(:) = 0.D0
         config_nematic_order_squared(:) = 0.D0
         config_stripe_x_binder(:) = 0.D0
         config_stripe_x_binder_squared(:) = 0.D0
         config_stripe_y_binder(:) = 0.D0
         config_stripe_y_binder_squared(:) = 0.D0
         config_local_nematic_mean(:) = 0.D0
         config_abs_local_nematic_mean(:) = 0.D0
         config_local_nematic_mean_squared(:) = 0.D0
         config_local_nematic_second_moment(:) = 0.D0
         config_local_nematic_fourth_moment(:) = 0.D0
         config_local_nematic_binder(:) = 0.D0
         config_local_nematic_binder_squared(:) = 0.D0

         config_energy_histogram(:, :) = 0.D0
         config_phi(:) = 0.D0
         config_pt_swap_attempts(:) = 0
         config_pt_swap_accepts(:) = 0
         config_pt_visits_both_ends = 0
         spin_sum_by_site(:) = 0.D0
         local_nematic_sum_by_site(:) = 0.D0
         local_nematic_sample_count = 0

      end subroutine reset_current_config_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine validate_run_options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (SAMPLERMODE /= 'LOCAL' .and. SAMPLERMODE /= 'PT') then
            print *, 'ERROR: SAMPLERMODE must be LOCAL or PT. Got ', SAMPLERMODE
            stop 1
         end if

         if (UPDATEORDER /= 'SEQ' .and. UPDATEORDER /= 'RND') then
            print *, 'ERROR: UPDATEORDER must be SEQ or RND. Got ', UPDATEORDER
            stop 1
         end if

         if (TEMPPROTOCOL /= 'CARRY' .and. TEMPPROTOCOL /= 'INDEP') then
            print *, 'ERROR: TEMPPROTOCOL must be CARRY or INDEP. Got ', TEMPPROTOCOL
            stop 1
         end if

         if (PTSTARTMODE /= 'COMMON' .and. PTSTARTMODE /= 'INDEP') then
            print *, 'ERROR: PTSTARTMODE must be COMMON or INDEP. Got ', PTSTARTMODE
            stop 1
         end if

         if (PTLADDERMODE /= 'LINEAR' .and. PTLADDERMODE /= 'SPLIT') then
            print *, 'ERROR: PTLADDERMODE must be LINEAR or SPLIT. Got ', PTLADDERMODE
            stop 1
         end if

         if (TMAX < TMIN) then
            print *, 'ERROR: TMAX must be >= TMIN. Got ', TMAX, TMIN
            stop 1
         end if

         if (SAMPLERMODE /= 'PT') then
            if (DT <= 0.D0) then
               print *, 'ERROR: DT must be positive. Got ', DT
               stop 1
            end if
         end if

         if (PTPARALLELMODE /= 'SERIAL' .and. PTPARALLELMODE /= 'GROUPED') then
            print *, 'ERROR: PTPARALLELMODE must be SERIAL or GROUPED. Got ', PTPARALLELMODE
            stop 1
         end if

         if (SAMPLERMODE == 'PT' .and. PTLADDERMODE == 'SPLIT') then
            if (NTEMP < 3) then
               print *, 'ERROR: SPLIT ladder requires at least 3 temperatures. Got ', NTEMP
               stop 1
            end if
            if (TLADDERMID <= TMIN .or. TLADDERMID >= TMAX) then
               print *, 'ERROR: TLADDERMID must lie strictly inside (TMIN,TMAX). Got ', TLADDERMID
               stop 1
            end if
            if (TLADDERFRAC <= 0.D0 .or. TLADDERFRAC >= 1.D0) then
               print *, 'ERROR: TLADDERFRAC must lie strictly inside (0,1). Got ', TLADDERFRAC
               stop 1
            end if
         end if

         if (SAMPLERMODE == 'PT') then
            if (NREPLICA /= NTEMP) then
               print *, 'ERROR: PT v1 requires NREPLICA = NTEMP. Got ', NREPLICA, NTEMP
               stop 1
            end if
            if (PTSWAPEVER <= 0) then
               print *, 'ERROR: PTSWAPEVER must be positive. Got ', PTSWAPEVER
               stop 1
            end if
            if (PTMEASINDEX < 1 .or. PTMEASINDEX > NTEMP) then
               print *, 'ERROR: PTMEASINDEX out of range. Got ', PTMEASINDEX
               stop 1
            end if
            if (PTPARALLELMODE == 'GROUPED') then
               if (WRITEHIST .eqv. .true. .or. WRITENUMFLIPS .eqv. .true. .or. WRITEBINDTIME .eqv. .true.) then
                  print *, 'ERROR: GROUPED PT v1 requires WRITEHIST, WRITENUMFLIPS, and WRITEBINDTIME to be false.'
                  stop 1
               end if
            end if
         else
            if (PTPARALLELMODE == 'GROUPED') then
               print *, 'ERROR: PTPARALLELMODE=GROUPED requires SAMPLERMODE=PT.'
               stop 1
            end if
         end if
      end subroutine validate_run_options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_temperature_ladder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: interval_index, ladder_index
         integer(i4b) :: interval_count, low_interval_count, high_interval_count

         interval_count = NTEMP - 1
         if (SAMPLERMODE == 'PT' .and. PTLADDERMODE == 'SPLIT' .and. interval_count > 0) then
            low_interval_count = nint(interval_count*TLADDERFRAC)
            low_interval_count = max(1, min(interval_count - 1, low_interval_count))
            high_interval_count = interval_count - low_interval_count

            temp_ladder(1) = TMIN
            do interval_index = 1, low_interval_count
               temp_ladder(interval_index + 1) = TMIN + &
     &            (TLADDERMID - TMIN)*real(interval_index, r8b)/real(low_interval_count, r8b)
            end do
            do interval_index = 1, high_interval_count
               temp_ladder(low_interval_count + 1 + interval_index) = TLADDERMID + &
     &            (TMAX - TLADDERMID)*real(interval_index, r8b)/real(high_interval_count, r8b)
            end do
         elseif (SAMPLERMODE == 'PT') then
            if (NTEMP <= 1) then
               temp_ladder(1) = TMIN
            else
               do ladder_index = 1, NTEMP
                  temp_ladder(ladder_index) = TMIN + &
     &               (TMAX - TMIN)*real(ladder_index - 1, r8b)/real(NTEMP - 1, r8b)
               end do
            end if
         else
            do ladder_index = 1, NTEMP
               temp_ladder(ladder_index) = TMIN + (ladder_index - 1)*DT
            end do
         end if

         do ladder_index = 1, NTEMP
            beta_ladder(ladder_index) = 1.D0/temp_ladder(ladder_index)
         end do
      end subroutine setup_temperature_ladder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_stripe_parities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: site2d, ix, iy

         do isite = 0, L3 - 1
            site2d = mod(isite, L2)
            ix = mod(site2d, L)
            iy = site2d/L
            if (mod(ix, 2) == 0) then
               stripe_parity_x(isite) = 1.D0
            else
               stripe_parity_x(isite) = -1.D0
            end if
            if (mod(iy, 2) == 0) then
               stripe_parity_y(isite) = 1.D0
            else
               stripe_parity_y(isite) = -1.D0
            end if
         end do
      end subroutine setup_stripe_parities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function temperature_for_index(temperature_index)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: temperature_index

         if (SAMPLERMODE == 'PT') then
            temperature_for_index = temp_ladder(temperature_index)
         elseif (START == 'HOT') then
            temperature_for_index = temp_ladder(NTEMP - temperature_index + 1)
         else
            temperature_for_index = temp_ladder(temperature_index)
         end if
      end function temperature_for_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reset_local_nematic_average
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         spin_sum_by_site(:) = 0.D0
         local_nematic_sum_by_site(:) = 0.D0
         local_nematic_sample_count = 0
      end subroutine reset_local_nematic_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function local_nematic_at_site(spins, site)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         integer(i4b), intent(in) :: site

         local_nematic_at_site = real(spins(site)*spins(m3(site))*occu(m3(site)), r8b) - &
                                 real(spins(site)*spins(m1(site))*occu(m1(site)), r8b) + &
                                 real(spins(site)*spins(m4(site))*occu(m4(site)), r8b) - &
                                 real(spins(site)*spins(m2(site))*occu(m2(site)), r8b)
      end function local_nematic_at_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine accumulate_local_nematic_average(spins)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)

         do isite = 0, L3 - 1
            if (occu(isite) == 1) then
               spin_sum_by_site(isite) = spin_sum_by_site(isite) + real(spins(isite), r8b)
               local_nematic_sum_by_site(isite) = local_nematic_sum_by_site(isite) + &
     &            local_nematic_at_site(spins, isite)
            end if
         end do
         local_nematic_sample_count = local_nematic_sample_count + 1
      end subroutine accumulate_local_nematic_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine spatial_mean_local_nematic(spins, mean_local_nematic)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         real(r8b), intent(out)   :: mean_local_nematic
         integer(i4b)             :: occ_count
         real(r8b)                :: local_nematic_sum

         local_nematic_sum = 0.D0
         occ_count = max(1, L3 - nimp)

         do isite = 0, L3 - 1
            if (occu(isite) == 0) cycle
            local_nematic_sum = local_nematic_sum + local_nematic_at_site(spins, isite)
         end do

         mean_local_nematic = local_nematic_sum/real(occ_count, r8b)
      end subroutine spatial_mean_local_nematic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reset_temperature_measurement_sums(sums)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         type(temperature_measurement_sums), intent(inout) :: sums

         sums%abs_magnetization(:) = 0.D0
         sums%magnetization_second_moment(:) = 0.D0
         sums%magnetization_fourth_moment(:) = 0.D0
         sums%energy(:) = 0.D0
         sums%energy_second_moment(:) = 0.D0
         sums%abs_stripe_x(:) = 0.D0
         sums%abs_stripe_y(:) = 0.D0
         sums%stripe_x_second_moment(:) = 0.D0
         sums%stripe_x_fourth_moment(:) = 0.D0
         sums%stripe_y_second_moment(:) = 0.D0
         sums%stripe_y_fourth_moment(:) = 0.D0
         sums%stripe_second_moment(:) = 0.D0
         sums%stripe_fourth_moment(:) = 0.D0
         sums%abs_stripe_order(:) = 0.D0
         sums%nematic_order(:) = 0.D0
         sums%local_nematic(:) = 0.D0
         sums%abs_local_nematic(:) = 0.D0
         sums%local_nematic_second_moment(:) = 0.D0
         sums%local_nematic_fourth_moment(:) = 0.D0
         sums%best_energy(:) = huge(1.D0)
         sums%final_energy(:) = 0.D0
      end subroutine reset_temperature_measurement_sums

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine accumulate_temperature_measurement(temperature_index, observables, mean_local_nematic, sums)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: temperature_index
         type(observable_state), intent(in) :: observables
         real(r8b), intent(in)    :: mean_local_nematic
         type(temperature_measurement_sums), intent(inout) :: sums
         real(r8b)                :: stripe_second_moment, mean_abs_stripe_order, mean_nematic_order

         ! Add one measured state into the bucket for its current temperature.
         stripe_second_moment = observables%stripe_x**2 + observables%stripe_y**2
         mean_abs_stripe_order = sqrt(stripe_second_moment)
         mean_nematic_order = abs(observables%stripe_y**2 - observables%stripe_x**2)

         sums%abs_magnetization(temperature_index) = sums%abs_magnetization(temperature_index) + abs(observables%magnetization)
         sums%magnetization_second_moment(temperature_index) = sums%magnetization_second_moment(temperature_index) + observables%magnetization**2
         sums%magnetization_fourth_moment(temperature_index) = sums%magnetization_fourth_moment(temperature_index) + observables%magnetization**4
         sums%energy(temperature_index) = sums%energy(temperature_index) + observables%energy
         sums%energy_second_moment(temperature_index) = sums%energy_second_moment(temperature_index) + observables%energy**2
         sums%abs_stripe_x(temperature_index) = sums%abs_stripe_x(temperature_index) + abs(observables%stripe_x)
         sums%abs_stripe_y(temperature_index) = sums%abs_stripe_y(temperature_index) + abs(observables%stripe_y)
         sums%stripe_x_second_moment(temperature_index) = sums%stripe_x_second_moment(temperature_index) + observables%stripe_x**2
         sums%stripe_x_fourth_moment(temperature_index) = sums%stripe_x_fourth_moment(temperature_index) + observables%stripe_x**4
         sums%stripe_y_second_moment(temperature_index) = sums%stripe_y_second_moment(temperature_index) + observables%stripe_y**2
         sums%stripe_y_fourth_moment(temperature_index) = sums%stripe_y_fourth_moment(temperature_index) + observables%stripe_y**4
         sums%stripe_second_moment(temperature_index) = sums%stripe_second_moment(temperature_index) + stripe_second_moment
         sums%stripe_fourth_moment(temperature_index) = sums%stripe_fourth_moment(temperature_index) + stripe_second_moment**2
         sums%abs_stripe_order(temperature_index) = sums%abs_stripe_order(temperature_index) + mean_abs_stripe_order
         sums%nematic_order(temperature_index) = sums%nematic_order(temperature_index) + mean_nematic_order
         sums%local_nematic(temperature_index) = sums%local_nematic(temperature_index) + mean_local_nematic
         sums%abs_local_nematic(temperature_index) = sums%abs_local_nematic(temperature_index) + abs(mean_local_nematic)
         sums%local_nematic_second_moment(temperature_index) = sums%local_nematic_second_moment(temperature_index) + mean_local_nematic**2
         sums%local_nematic_fourth_moment(temperature_index) = sums%local_nematic_fourth_moment(temperature_index) + mean_local_nematic**4
         sums%best_energy(temperature_index) = min(sums%best_energy(temperature_index), observables%energy)
      end subroutine accumulate_temperature_measurement

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine summarize_one_measured_state(observables, mean_local_nematic, summary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         type(observable_state), intent(in) :: observables
         real(r8b), intent(in) :: mean_local_nematic
         type(config_observable_summary), intent(out) :: summary
         real(r8b) :: stripe_second_moment

         stripe_second_moment = observables%stripe_x**2 + observables%stripe_y**2

         summary%mean_abs_magnetization = abs(observables%magnetization)
         summary%magnetization_second_moment = observables%magnetization**2
         summary%magnetization_fourth_moment = observables%magnetization**4
         summary%mean_energy = observables%energy
         summary%energy_second_moment = observables%energy**2
         summary%mean_abs_stripe_x = abs(observables%stripe_x)
         summary%mean_abs_stripe_y = abs(observables%stripe_y)
         summary%stripe_x_second_moment = observables%stripe_x**2
         summary%stripe_x_fourth_moment = observables%stripe_x**4
         summary%stripe_y_second_moment = observables%stripe_y**2
         summary%stripe_y_fourth_moment = observables%stripe_y**4
         summary%stripe_second_moment = stripe_second_moment
         summary%stripe_fourth_moment = stripe_second_moment**2
         summary%mean_abs_stripe_order = sqrt(stripe_second_moment)
         summary%mean_nematic_order = abs(observables%stripe_y**2 - observables%stripe_x**2)
         summary%mean_local_nematic = mean_local_nematic
         summary%mean_abs_local_nematic = abs(mean_local_nematic)
         summary%local_nematic_second_moment = mean_local_nematic**2
         summary%local_nematic_fourth_moment = mean_local_nematic**4
         summary%best_energy = observables%energy
         summary%final_energy = observables%energy
      end subroutine summarize_one_measured_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine average_measurement_sums_for_temperature(sums, temperature_index, nmeasure, summary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         type(temperature_measurement_sums), intent(in) :: sums
         integer(i4b), intent(in) :: temperature_index, nmeasure
         type(config_observable_summary), intent(out) :: summary
         real(r8b) :: measurement_count

         measurement_count = real(nmeasure, r8b)

         summary%mean_abs_magnetization = sums%abs_magnetization(temperature_index)/measurement_count
         summary%magnetization_second_moment = sums%magnetization_second_moment(temperature_index)/measurement_count
         summary%magnetization_fourth_moment = sums%magnetization_fourth_moment(temperature_index)/measurement_count
         summary%mean_energy = sums%energy(temperature_index)/measurement_count
         summary%energy_second_moment = sums%energy_second_moment(temperature_index)/measurement_count
         summary%mean_abs_stripe_x = sums%abs_stripe_x(temperature_index)/measurement_count
         summary%mean_abs_stripe_y = sums%abs_stripe_y(temperature_index)/measurement_count
         summary%stripe_x_second_moment = sums%stripe_x_second_moment(temperature_index)/measurement_count
         summary%stripe_x_fourth_moment = sums%stripe_x_fourth_moment(temperature_index)/measurement_count
         summary%stripe_y_second_moment = sums%stripe_y_second_moment(temperature_index)/measurement_count
         summary%stripe_y_fourth_moment = sums%stripe_y_fourth_moment(temperature_index)/measurement_count
         summary%stripe_second_moment = sums%stripe_second_moment(temperature_index)/measurement_count
         summary%stripe_fourth_moment = sums%stripe_fourth_moment(temperature_index)/measurement_count
         summary%mean_abs_stripe_order = sums%abs_stripe_order(temperature_index)/measurement_count
         summary%mean_nematic_order = sums%nematic_order(temperature_index)/measurement_count
         summary%mean_local_nematic = sums%local_nematic(temperature_index)/measurement_count
         summary%mean_abs_local_nematic = sums%abs_local_nematic(temperature_index)/measurement_count
         summary%local_nematic_second_moment = sums%local_nematic_second_moment(temperature_index)/measurement_count
         summary%local_nematic_fourth_moment = sums%local_nematic_fourth_moment(temperature_index)/measurement_count
         summary%best_energy = sums%best_energy(temperature_index)
         summary%final_energy = sums%final_energy(temperature_index)
      end subroutine average_measurement_sums_for_temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function stripe_coexistence_angle(stripe_x, stripe_y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b), intent(in) :: stripe_x, stripe_y
         real(r8b)             :: denom

         denom = stripe_x**2 + stripe_y**2
         if (denom > 0.D0) then
            ! This equals 0.5*(1-cos(4 theta)): 1 on the diagonal and 0 on either axis.
            stripe_coexistence_angle = 4.D0*stripe_x**2*stripe_y**2/(denom**2)
         else
            stripe_coexistence_angle = 0.D0
         end if
      end function stripe_coexistence_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine validate_observable_state(spins, observables, context, diagnostic_index_1, diagnostic_index_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         type(observable_state), intent(in) :: observables
         character(len=*), intent(in) :: context
         integer(i4b), intent(in) :: diagnostic_index_1, diagnostic_index_2
         type(observable_state)   :: reference_observables
         real(r8b), parameter     :: tol = 1.0D-10

         if (.not. VALIDATE_INCREMENTAL) return

         call measure_observables(spins, reference_observables)
         if (abs(observables%magnetization - reference_observables%magnetization) > tol .or. &
     &       abs(observables%energy - reference_observables%energy) > tol .or. &
     &       abs(observables%stripe_x - reference_observables%stripe_x) > tol .or. &
     &       abs(observables%stripe_y - reference_observables%stripe_y) > tol) then
            print *, 'ERROR: incremental observable mismatch in ', trim(context), diagnostic_index_1, diagnostic_index_2
            print *, ' current=', observables%magnetization, observables%energy, observables%stripe_x, observables%stripe_y
            print *, ' reference=', reference_observables%magnetization, reference_observables%energy, &
     &         reference_observables%stripe_x, reference_observables%stripe_y
            stop 1
         end if
      end subroutine validate_observable_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine prepare_temperature_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: temp_seed

         if (TEMPPROTOCOL == 'INDEP') then
            ! Reset the RNG before each temperature so fresh-start tests are reproducible
            ! and independent of how many random numbers earlier temperatures consumed.
            temp_seed = IRINIT + (iconf - 1)*NTEMP + temperature_index
            call kissinit(temp_seed)
            call initialize_spin_array(s1)
         end if
      end subroutine prepare_temperature_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_vacancies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         integer(i4b) :: trials, max_trials

         if (IMPMODE == 'RAN') then                ! uncorrelated impurities
            nimp = 0
            do isite = 0, L3 - 1
               if (rkiss05() < PP) then
                  occu(isite) = 0
                  nimp = nimp + 1
               else
                  occu(isite) = 1
               end if
            end do
         end if

         if (IMPMODE == 'COR') then                ! anticorrelated impurities
            occu(:) = 1
            nimp = 0
            max_trials = 25*L3
            trials = 0
            do while (nimp < PP*L3 .and. trials < max_trials)
               isite = int(L3*rkiss05())
               if (occu(isite)*occu(m1(isite))*occu(m2(isite))*occu(m3(isite))*occu(m4(isite)) == 1) then
                  occu(isite) = 0
                  nimp = nimp + 1
               end if
               trials = trials + 1
            end do
            if (nimp < PP*L3) then
               print *, "Warning: could not add all impurities. Only added", nimp
            end if
         end if

         do isite = 0, L3 - 1
            if (rkiss05() < PPh) then
               jh(isite) = 0
            else
               jh(isite) = 1
            end if
         end do

         do isite = 0, L3 - 1
            if (rkiss05() < PPv) then
               jv(isite) = 0
            else
               jv(isite) = 1
            end if
         end do

      end subroutine setup_vacancies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mark_vacancy_flipcounts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: temp_slot

         do temp_slot = 1, NTEMP
            where (occu(:) == 0)
               flip_counts_this_config(temp_slot, :) = -1
            end where
         end do
      end subroutine mark_vacancy_flipcounts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialize_spin_array(spins)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(out) :: spins(0:L3 - 1)
         integer(i4b)    :: ix, iy, iz

         if (START == 'HOT') then
            do isite = 0, L3 - 1
               if (rkiss05() < 0.5D0) then
                  spins(isite) = 1
               else
                  spins(isite) = -1
               end if
            end do
         elseif (START == 'FER') then
            spins(:) = 1
         else                                  !start from stripe configuration
            do iz = 0, H - 1
               do iy = 0, L - 1
                  do ix = 0, L - 1
                     isite = L*iy + ix + L2*iz
                     if (STRAXIS == 'X') then
                        spins(isite) = (-1)**iy
                     else
                        spins(isite) = (-1)**ix
                     end if
                  end do
               end do
            end do
         end if

      end subroutine initialize_spin_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialize_pt_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: replica_index
         integer(i4b) :: replica_seed

         ! Reset PT diagnostics for the new disorder sample.
         config_pt_swap_attempts(:) = 0
         config_pt_swap_accepts(:) = 0
         config_pt_visits_both_ends = 0
         config_energy_histogram(:, :) = 0.D0
         energy_histogram_min(:) = 0.D0
         energy_histogram_max(:) = 0.D0
         energy_histogram_bin_size(:) = 0.D0
         flip_counts_this_config(:, :) = 0
         call mark_vacancy_flipcounts

         ! Start with replica k sitting at temperature index k.
         do replica_index = 1, NREPLICA
            replica_at_temperature(replica_index) = replica_index
            temperature_of_replica(replica_index) = replica_index
            replica_seed = IRINIT + (iconf - 1)*NREPLICA + replica_index
            call kissinit(replica_seed)
            ! Each replica gets its own spin configuration on the same disorder sample.
            call initialize_spin_array(s1_rep(0:L3 - 1, replica_index))
         end do

         ! Use a dedicated post-initialization seed so PT dynamics remain
         ! reproducible even when the number of replicas changes later.
         call kissinit(IRINIT + 100000 + iconf)
      end subroutine initialize_pt_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pt_equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep
         integer(i4b) :: swap_phase
         integer(i4b) :: replica_index
         integer(i4b) :: progress_every
         integer(kind=8) :: phase_clock_start
         type(observable_state) :: replica_observables(NREPLICA)

         call timer_start(phase_clock_start)
         do replica_index = 1, NREPLICA
            call measure_observables(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index))
            call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index), &
     &         'PTEQINIT', 0, replica_index)
         end do

         ! Alternate odd and even neighbor pairs on successive swap passes.
         swap_phase = 1
         progress_every = max(1, NEQ/10) ! write progress every 10% of the equilibration sweeps
         call write_progress('PT_EQUIL', 0, NEQ)
         do isweep = 1, NEQ
            call pt_cycle(swap_phase, replica_observables)
            ! Swap phases alternate between ladder bonds starting at 1 and 2.
            swap_phase = 3 - swap_phase
            if (mod(isweep, progress_every) == 0 .or. isweep == NEQ) call write_progress('PT_EQUIL', isweep, NEQ)
         end do
         call write_timing('PT_EQUIL', phase_clock_start)
      end subroutine pt_equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pt_cycle(swap_phase, replica_observables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: swap_phase
         type(observable_state), intent(inout) :: replica_observables(NREPLICA)
         integer(i4b)             :: temperature_index, replica_index, isweep

         ! Evolve each temperature index with whichever replica currently sits there.
         do temperature_index = 1, NREPLICA
            replica_index = replica_at_temperature(temperature_index)
            do isweep = 1, PTSWAPEVER
               call metropolis_sweep(s1_rep(0:L3 - 1, replica_index), beta_ladder(temperature_index), &
     &            flip_counts_this_config(temperature_index, :), replica_observables(replica_index))
            end do
            call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index), &
     &         'PTCYCLE', temperature_index, replica_index)
         end do

         ! Then try the active set of nearest-neighbor replica exchanges.
         call pt_attempt_swaps(swap_phase, replica_observables)
      end subroutine pt_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pt_attempt_swaps(swap_phase, replica_observables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: swap_phase
         type(observable_state), intent(in) :: replica_observables(NREPLICA)
         integer(i4b)             :: temperature_index, swap_start_index
         integer(i4b)             :: lower_replica_index, upper_replica_index, saved_lower_replica_index
         real(r8b)                :: delta
         real(r8b)                :: lower_total_energy, upper_total_energy

         if (swap_phase == 1) then
            swap_start_index = 1
         else
            swap_start_index = 2
         end if

         ! Swap neighboring ladder slots only, alternating odd/even bonds.
         do temperature_index = swap_start_index, NREPLICA - 1, 2
            lower_replica_index = replica_at_temperature(temperature_index)
            upper_replica_index = replica_at_temperature(temperature_index + 1)
            config_pt_swap_attempts(temperature_index) = config_pt_swap_attempts(temperature_index) + 1

            lower_total_energy = replica_observables(lower_replica_index)%energy*(L3 - nimp)
            upper_total_energy = replica_observables(upper_replica_index)%energy*(L3 - nimp)
            delta = (beta_ladder(temperature_index) - beta_ladder(temperature_index + 1))*(lower_total_energy - upper_total_energy)

            ! Accept with the standard replica-exchange Metropolis rule.
            if (delta >= 0.D0 .or. rkiss05() < exp(delta)) then
               config_pt_swap_accepts(temperature_index) = config_pt_swap_accepts(temperature_index) + 1
               saved_lower_replica_index = replica_at_temperature(temperature_index)
               replica_at_temperature(temperature_index) = replica_at_temperature(temperature_index + 1)
               replica_at_temperature(temperature_index + 1) = saved_lower_replica_index
               ! Keep the inverse map in sync after the swap.
               temperature_of_replica(lower_replica_index) = temperature_index + 1
               temperature_of_replica(upper_replica_index) = temperature_index
            end if
         end do
      end subroutine pt_attempt_swaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pt_measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, temperature_index, replica_index, nmeasure
         integer(i4b) :: energy_histogram_bin_index
         integer(i4b) :: swap_phase
         integer(i4b) :: progress_every
         integer(kind=8) :: phase_clock_start
         type(observable_state) :: replica_observables(NREPLICA)
         type(config_observable_summary) :: summary
         type(temperature_measurement_sums) :: measurement_sums
         real(r8b)    :: mean_local_nematic_by_replica(NREPLICA)
         real(r8b)    :: measured_energy, measured_local_nematic
         real(r8b)    :: target_stripe_x_second_moment_sum, target_stripe_x_fourth_moment_sum
         real(r8b)    :: target_stripe_y_second_moment_sum, target_stripe_y_fourth_moment_sum
         real(r8b)    :: mean_stripe_x_second_moment, mean_stripe_x_fourth_moment
         real(r8b)    :: mean_stripe_y_second_moment, mean_stripe_y_fourth_moment
         real(r8b)    :: stripe_x_binder_value, stripe_y_binder_value
         logical(ilog) :: do_bind_time, visited_cold(NREPLICA), visited_hot(NREPLICA), do_trace
         character(len=21) :: bt_trace_file, trace_file

         call reset_temperature_measurement_sums(measurement_sums)
         config_phi(:) = p(:)
         config_flip_counts(:, :) = flip_counts_this_config(:, :)
         call reset_local_nematic_average
         call timer_start(phase_clock_start)

         ! Track whether a replica has visited both ends of the ladder.
         visited_cold(:) = .false.
         visited_hot(:) = .false.
         do_bind_time = WRITEBINDTIME .and. (iconf <= NTRACECONF)  ! only trace a few configs so files do not explode
         do_trace = (iconf <= NTRACECONF)

         if (do_bind_time) then
            ! This is a running Binder trace at the selected PT temperature.
            write (bt_trace_file, '("cfbt",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (77, file=bt_trace_file, status='replace')
            write (77, *) 'iconf T'
            write (77, '(I4,1x,E12.6)') iconf, temp_ladder(PTMEASINDEX)
            write (77, *) 'isweep, binderx, bindery'
            target_stripe_x_second_moment_sum = 0.D0
            target_stripe_x_fourth_moment_sum = 0.D0
            target_stripe_y_second_moment_sum = 0.D0
            target_stripe_y_fourth_moment_sum = 0.D0
         end if

         if (do_trace) then
            ! Track a few replica positions to make sure swaps are actually moving.
            write (trace_file, '("cftr",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (78, file=trace_file, status='replace')
            write (78, *) 'isweep replica temp_index T'
         end if

         if (NMESS <= 0) then
            call write_progress('PT_MEAS', 0, 0)
            do replica_index = 1, NREPLICA
               ! With zero measurement sweeps, just measure the current replica state once.
               call measure_observables(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index))
               call spatial_mean_local_nematic(s1_rep(0:L3 - 1, replica_index), mean_local_nematic_by_replica(replica_index))
               call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index), &
     &            'PTINIT', 0, replica_index)
            end do
            call accumulate_local_nematic_average(s1_rep(0:L3 - 1, replica_at_temperature(PTMEASINDEX)))
            nmeasure = 1
         else
            ! In PT mode NEQ and NMESS count whole PT cycles, not raw sweeps.
            nmeasure = NMESS
            swap_phase = 1
            progress_every = max(1, NMESS/10)
            do replica_index = 1, NREPLICA
               call measure_observables(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index))
               call spatial_mean_local_nematic(s1_rep(0:L3 - 1, replica_index), mean_local_nematic_by_replica(replica_index))
               call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables(replica_index), &
     &            'PTINIT', 0, replica_index)
            end do
            call write_progress('PT_MEAS', 0, NMESS)
            do isweep = 1, NMESS
               call pt_cycle(swap_phase, replica_observables)  ! local sweeps, then one swap parity

               do replica_index = 1, NREPLICA
                  ! The spin array follows the replica, not the temperature index.
                  if (temperature_of_replica(replica_index) == 1) visited_cold(replica_index) = .true.
                  if (temperature_of_replica(replica_index) == NREPLICA) visited_hot(replica_index) = .true.
                  call spatial_mean_local_nematic(s1_rep(0:L3 - 1, replica_index), mean_local_nematic_by_replica(replica_index))
               end do

               if (do_trace) then
                  do replica_index = 1, min(NREPLICA, 3)
                     write (78, '(I8,1x,I4,1x,I4,1x,E12.6)') isweep, replica_index, temperature_of_replica(replica_index), temp_ladder(temperature_of_replica(replica_index))
                  end do
               end if

               if (do_bind_time) then
                  replica_index = replica_at_temperature(PTMEASINDEX)
                  target_stripe_x_second_moment_sum = target_stripe_x_second_moment_sum + replica_observables(replica_index)%stripe_x**2
                  target_stripe_x_fourth_moment_sum = target_stripe_x_fourth_moment_sum + replica_observables(replica_index)%stripe_x**4
                  target_stripe_y_second_moment_sum = target_stripe_y_second_moment_sum + replica_observables(replica_index)%stripe_y**2
                  target_stripe_y_fourth_moment_sum = target_stripe_y_fourth_moment_sum + replica_observables(replica_index)%stripe_y**4

                  mean_stripe_x_second_moment = target_stripe_x_second_moment_sum/isweep
                  mean_stripe_x_fourth_moment = target_stripe_x_fourth_moment_sum/isweep
                  mean_stripe_y_second_moment = target_stripe_y_second_moment_sum/isweep
                  mean_stripe_y_fourth_moment = target_stripe_y_fourth_moment_sum/isweep

                  if (mean_stripe_x_second_moment > 0.D0) then
                     stripe_x_binder_value = 1.D0 - &
     &                  mean_stripe_x_fourth_moment/(3.D0*mean_stripe_x_second_moment**2)
                  else
                     stripe_x_binder_value = 0.D0
                  end if
                  if (mean_stripe_y_second_moment > 0.D0) then
                     stripe_y_binder_value = 1.D0 - &
     &                  mean_stripe_y_fourth_moment/(3.D0*mean_stripe_y_second_moment**2)
                  else
                     stripe_y_binder_value = 0.D0
                  end if
                  write (77, '(I8,1x,E12.6,1x,E12.6)') isweep, stripe_x_binder_value, stripe_y_binder_value
               end if

               call accumulate_local_nematic_average(s1_rep(0:L3 - 1, replica_at_temperature(PTMEASINDEX)))

               ! Accumulate observables by temperature index after swaps.
               do temperature_index = 1, NTEMP
                  replica_index = replica_at_temperature(temperature_index)             ! after swaps, find which replica is at this T
                  measured_energy = replica_observables(replica_index)%energy
                  measured_local_nematic = mean_local_nematic_by_replica(replica_index)
                  call accumulate_temperature_measurement(temperature_index, replica_observables(replica_index), &
     &               measured_local_nematic, measurement_sums)
                  measurement_sums%final_energy(temperature_index) = measured_energy

                  if (WRITEHIST .eqv. .true.) then
                     if (energy_histogram_bin_size(temperature_index) <= 0.D0) then
                        energy_histogram_min(temperature_index) = measured_energy - (abs(measured_energy) + 1.D0)
                        energy_histogram_max(temperature_index) = measured_energy + (abs(measured_energy) + 1.D0)
                        energy_histogram_bin_size(temperature_index) = (energy_histogram_max(temperature_index) - energy_histogram_min(temperature_index))/BINS
                     end if
                     energy_histogram_bin_index = 1 + int((measured_energy - energy_histogram_min(temperature_index)) &
     &                  /energy_histogram_bin_size(temperature_index))
                     energy_histogram_bin_index = max(1, min(BINS, energy_histogram_bin_index))
                     config_energy_histogram(temperature_index, energy_histogram_bin_index) = &
     &                  config_energy_histogram(temperature_index, energy_histogram_bin_index) + 1.D0
                  end if
               end do

               swap_phase = 3 - swap_phase
               if (mod(isweep, progress_every) == 0 .or. isweep == NMESS) call write_progress('PT_MEAS', isweep, NMESS)
            end do
         end if
         call write_timing('PT_MEAS', phase_clock_start)

         if (do_bind_time) close (77)
         if (do_trace) close (78)

         ! Count replicas that completed at least one cold-to-hot visitation.
         config_pt_visits_both_ends = count(visited_cold(:) .and. visited_hot(:))

         do temperature_index = 1, NTEMP
            replica_index = replica_at_temperature(temperature_index)
            if (NMESS <= 0) then
               ! No sweep list exists here, so use the one measured value directly.
               call summarize_one_measured_state(replica_observables(replica_index), mean_local_nematic_by_replica(replica_index), summary)
            else
               ! Convert the accumulated sums into thermal averages for this config.
               call average_measurement_sums_for_temperature(measurement_sums, temperature_index, nmeasure, summary)
            end if

            call store_config_observables(temperature_index, beta_ladder(temperature_index), summary)

            if (WRITEHIST .eqv. .true. .and. NMESS > 0) then
               config_energy_histogram(temperature_index, :) = config_energy_histogram(temperature_index, :)/NMESS
            end if
            config_flip_counts(temperature_index, :) = flip_counts_this_config(temperature_index, :)
         end do
      end subroutine pt_measurement_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      subroutine grouped_pt_run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (.not. pt_group_active) return
         call grouped_initialize_pt_state
         call grouped_pt_equilibration_loop
         call grouped_pt_measurement_loop
      end subroutine grouped_pt_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reduce_grouped_pt_measurement_sums(local_sums, group_sums)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         type(temperature_measurement_sums), intent(in) :: local_sums
         type(temperature_measurement_sums), intent(out) :: group_sums

         call MPI_REDUCE(local_sums%abs_magnetization, group_sums%abs_magnetization, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%magnetization_second_moment, group_sums%magnetization_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%magnetization_fourth_moment, group_sums%magnetization_fourth_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%energy, group_sums%energy, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%energy_second_moment, group_sums%energy_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%abs_stripe_x, group_sums%abs_stripe_x, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%abs_stripe_y, group_sums%abs_stripe_y, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_x_second_moment, group_sums%stripe_x_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_x_fourth_moment, group_sums%stripe_x_fourth_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_y_second_moment, group_sums%stripe_y_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_y_fourth_moment, group_sums%stripe_y_fourth_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_second_moment, group_sums%stripe_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%stripe_fourth_moment, group_sums%stripe_fourth_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%abs_stripe_order, group_sums%abs_stripe_order, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%nematic_order, group_sums%nematic_order, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%local_nematic, group_sums%local_nematic, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%abs_local_nematic, group_sums%abs_local_nematic, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%local_nematic_second_moment, group_sums%local_nematic_second_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%local_nematic_fourth_moment, group_sums%local_nematic_fourth_moment, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%best_energy, group_sums%best_energy, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_MIN, 0, pt_comm, ierr)
         call MPI_REDUCE(local_sums%final_energy, group_sums%final_energy, &
     &      NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
      end subroutine reduce_grouped_pt_measurement_sums

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grouped_initialize_pt_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: replica_index
         integer(i4b) :: replica_seed

         config_pt_swap_attempts(:) = 0
         config_pt_swap_accepts(:) = 0
         config_pt_visits_both_ends = 0
         config_energy_histogram(:, :) = 0.D0
         energy_histogram_min(:) = 0.D0
         energy_histogram_max(:) = 0.D0
         energy_histogram_bin_size(:) = 0.D0
         flip_counts_this_config(:, :) = 0
         call mark_vacancy_flipcounts

         if (pt_group_master) then
            do replica_index = 1, NREPLICA
               replica_at_temperature(replica_index) = replica_index
               temperature_of_replica(replica_index) = replica_index
            end do
         end if
         call MPI_BCAST(replica_at_temperature, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)
         call MPI_BCAST(temperature_of_replica, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)

         replica_index = pt_group_rank + 1
         replica_seed = IRINIT + (iconf - 1)*NREPLICA + replica_index
         call kissinit(replica_seed)
         call initialize_spin_array(s1_rep(0:L3 - 1, replica_index))
         call kissinit(IRINIT + 100000 + (iconf - 1)*NREPLICA + replica_index)
      end subroutine grouped_initialize_pt_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grouped_pt_equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, swap_phase, replica_index, current_temperature_index, progress_every
         integer(kind=8) :: phase_clock_start
         type(observable_state) :: replica_observables

         call timer_start(phase_clock_start)
         replica_index = pt_group_rank + 1
         current_temperature_index = temperature_of_replica(replica_index)
         call measure_observables(s1_rep(0:L3 - 1, replica_index), replica_observables)
         call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables, &
     &      'GPTINIT', current_temperature_index, replica_index)

         swap_phase = 1
         progress_every = max(1, NEQ/10)
         if (pt_group_master) call write_progress('PT_EQUIL', 0, NEQ)
         do isweep = 1, NEQ
            call grouped_pt_cycle(swap_phase, replica_index, current_temperature_index, replica_observables)
            swap_phase = 3 - swap_phase
            if (pt_group_master .and. (mod(isweep, progress_every) == 0 .or. isweep == NEQ)) &
     &         call write_progress('PT_EQUIL', isweep, NEQ)
         end do
         call write_timing('PT_EQUIL', phase_clock_start)
      end subroutine grouped_pt_equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grouped_pt_cycle(swap_phase, replica_index, current_temperature_index, replica_observables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in)    :: swap_phase, replica_index
         integer(i4b), intent(inout) :: current_temperature_index
         type(observable_state), intent(inout) :: replica_observables
         integer(i4b)                :: isweep

         do isweep = 1, PTSWAPEVER
            call metropolis_sweep(s1_rep(0:L3 - 1, replica_index), beta_ladder(current_temperature_index), flip_counts_this_config(current_temperature_index, :), &
     &         replica_observables)
         end do
         call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables, &
     &      'GPTCYCLE', current_temperature_index, replica_index)
         call grouped_pt_attempt_swaps(swap_phase, replica_index, current_temperature_index, replica_observables%energy)
      end subroutine grouped_pt_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grouped_pt_attempt_swaps(swap_phase, replica_index, current_temperature_index, replica_energy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in)    :: swap_phase, replica_index
         integer(i4b), intent(inout) :: current_temperature_index
         real(r8b), intent(in)       :: replica_energy
         integer(i4b)                :: temperature_index, swap_start_index
         integer(i4b)                :: lower_replica_index, upper_replica_index, saved_lower_replica_index
         real(r8b)                   :: energy_by_replica(NREPLICA), delta, lower_total_energy, upper_total_energy

         energy_by_replica(:) = 0.D0
         call MPI_GATHER(replica_energy, 1, MPI_DOUBLE_PRECISION, energy_by_replica, 1, MPI_DOUBLE_PRECISION, 0, pt_comm, ierr)

         if (pt_group_master) then
            if (swap_phase == 1) then
               swap_start_index = 1
            else
               swap_start_index = 2
            end if
            do temperature_index = swap_start_index, NREPLICA - 1, 2
               lower_replica_index = replica_at_temperature(temperature_index)
               upper_replica_index = replica_at_temperature(temperature_index + 1)
               config_pt_swap_attempts(temperature_index) = config_pt_swap_attempts(temperature_index) + 1

               lower_total_energy = energy_by_replica(lower_replica_index)*(L3 - nimp)
               upper_total_energy = energy_by_replica(upper_replica_index)*(L3 - nimp)
               delta = (beta_ladder(temperature_index) - beta_ladder(temperature_index + 1))*(lower_total_energy - upper_total_energy)
               if (delta >= 0.D0 .or. rkiss05() < exp(delta)) then
                  config_pt_swap_accepts(temperature_index) = config_pt_swap_accepts(temperature_index) + 1
                  saved_lower_replica_index = replica_at_temperature(temperature_index)
                  replica_at_temperature(temperature_index) = replica_at_temperature(temperature_index + 1)
                  replica_at_temperature(temperature_index + 1) = saved_lower_replica_index
                  temperature_of_replica(lower_replica_index) = temperature_index + 1
                  temperature_of_replica(upper_replica_index) = temperature_index
               end if
            end do
         end if

         call MPI_BCAST(temperature_of_replica, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)
         current_temperature_index = temperature_of_replica(replica_index)
      end subroutine grouped_pt_attempt_swaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grouped_pt_measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, temperature_index, replica_index, current_temperature_index, nmeasure, swap_phase, progress_every
         integer(kind=8) :: phase_clock_start
         type(observable_state) :: replica_observables
         type(config_observable_summary) :: summary
         type(temperature_measurement_sums) :: local_measurement_sums, group_measurement_sums
         real(r8b)    :: mean_local_nematic_for_rank
         real(r8b)    :: local_visitboth(1), group_visitboth(1)
         real(r8b)    :: group_spin_sum_by_site(0:L3 - 1), group_local_nematic_sum_by_site(0:L3 - 1)
         real(r8b)    :: local_nematic_samples_for_rank(1), local_nematic_samples_for_group(1)
         logical(ilog) :: visited_cold, visited_hot

         call reset_temperature_measurement_sums(local_measurement_sums)
         config_phi(:) = p(:)
         call reset_local_nematic_average
         call timer_start(phase_clock_start)

         replica_index = pt_group_rank + 1                      ! in grouped PT each rank owns one replica
         current_temperature_index = temperature_of_replica(replica_index)
         visited_cold = .false.
         visited_hot = .false.
         call measure_observables(s1_rep(0:L3 - 1, replica_index), replica_observables)
         call spatial_mean_local_nematic(s1_rep(0:L3 - 1, replica_index), mean_local_nematic_for_rank)
         call validate_observable_state(s1_rep(0:L3 - 1, replica_index), replica_observables, &
     &      'GPTMEASINIT', current_temperature_index, replica_index)

         if (NMESS <= 0) then
            nmeasure = 1
            if (pt_group_master) call write_progress('PT_MEAS', 0, 0)
            if (WRITELOCNEM .eqv. .true. .and. current_temperature_index == PTMEASINDEX) &
     &         call accumulate_local_nematic_average(s1_rep(0:L3 - 1, replica_index))
            call accumulate_temperature_measurement(current_temperature_index, replica_observables, mean_local_nematic_for_rank, &
     &         local_measurement_sums)
         else
            nmeasure = NMESS
            swap_phase = 1
            progress_every = max(1, NMESS/10)
            if (pt_group_master) call write_progress('PT_MEAS', 0, NMESS)
            do isweep = 1, NMESS
               ! Grouped PT has to swap by talking inside pt_comm.
               call grouped_pt_cycle(swap_phase, replica_index, current_temperature_index, replica_observables)
               if (current_temperature_index == 1) visited_cold = .true.
               if (current_temperature_index == NREPLICA) visited_hot = .true.
               call spatial_mean_local_nematic(s1_rep(0:L3 - 1, replica_index), mean_local_nematic_for_rank)
               if (WRITELOCNEM .eqv. .true. .and. current_temperature_index == PTMEASINDEX) &
     &            call accumulate_local_nematic_average(s1_rep(0:L3 - 1, replica_index))
               call accumulate_temperature_measurement(current_temperature_index, replica_observables, mean_local_nematic_for_rank, &
     &            local_measurement_sums)

               swap_phase = 3 - swap_phase
               if (pt_group_master .and. (mod(isweep, progress_every) == 0 .or. isweep == NMESS)) &
     &            call write_progress('PT_MEAS', isweep, NMESS)
            end do
         end if
         call write_timing('PT_MEAS', phase_clock_start)

         local_measurement_sums%final_energy(:) = 0.D0
         local_measurement_sums%final_energy(current_temperature_index) = replica_observables%energy
         local_visitboth(1) = merge(1.D0, 0.D0, visited_cold .and. visited_hot)
         local_nematic_samples_for_rank(1) = real(local_nematic_sample_count, r8b)

         ! Put the per-replica sums back on the group master before storing them.
         call reduce_grouped_pt_measurement_sums(local_measurement_sums, group_measurement_sums)
         call MPI_REDUCE(local_visitboth, group_visitboth, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         if (WRITELOCNEM .eqv. .true.) then
            call MPI_REDUCE(spin_sum_by_site, group_spin_sum_by_site, L3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
            call MPI_REDUCE(local_nematic_sum_by_site, group_local_nematic_sum_by_site, L3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
            call MPI_REDUCE(local_nematic_samples_for_rank, local_nematic_samples_for_group, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         end if

         if (pt_group_master) then
            if (WRITELOCNEM .eqv. .true.) then
               ! Local nematic output is only written once, so gather its arrays here too.
               spin_sum_by_site(:) = group_spin_sum_by_site(:)
               local_nematic_sum_by_site(:) = group_local_nematic_sum_by_site(:)
               local_nematic_sample_count = nint(local_nematic_samples_for_group(1))
               temperature_index = PTMEASINDEX
               T = temp_ladder(PTMEASINDEX)
               beta = beta_ladder(PTMEASINDEX)
            end if
            config_pt_visits_both_ends = nint(group_visitboth(1))
            do temperature_index = 1, NTEMP
               ! The group sums are already by temperature index after the reductions.
               call average_measurement_sums_for_temperature(group_measurement_sums, temperature_index, nmeasure, summary)
               call store_config_observables(temperature_index, beta_ladder(temperature_index), summary)
            end do
         end if
      end subroutine grouped_pt_measurement_loop

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine snapshot_pt_target_replica
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: replica_index

         replica_index = replica_at_temperature(PTMEASINDEX)
         s1(:) = s1_rep(:, replica_index)
         temperature_index = PTMEASINDEX
         T = temp_ladder(PTMEASINDEX)
         beta = beta_ladder(PTMEASINDEX)
      end subroutine snapshot_pt_target_replica

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         type(observable_state) :: current_observables

         call measure_observables(s1, current_observables)
         call validate_observable_state(s1, current_observables, 'LOCALEQINIT', 0, temperature_index)

         do isweep = 1, NEQ
            call metropolis_sweep(s1, beta, flip_counts_this_config(temperature_index, :), current_observables)
            call validate_observable_state(s1, current_observables, 'LOCALEQ', isweep, temperature_index)
         end do

      end subroutine equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         type(observable_state) :: current_observables
         real(r8b)       :: current_local_nematic_mean
         real(r8b)       :: disorder_stripe_x_second_moment, disorder_stripe_x_fourth_moment, disorder_stripe_y_second_moment, disorder_stripe_y_fourth_moment
         real(r8b)       :: mean_stripe_x_second_moment, mean_stripe_x_fourth_moment
         real(r8b)       :: mean_stripe_y_second_moment, mean_stripe_y_fourth_moment
         real(r8b)       :: stripe_x_binder_value, stripe_y_binder_value
         logical(ilog)   :: do_bind_time
         character(len=21) :: bt_trace_file

         call reset_local_nematic_average
         call measure_observables(s1, current_observables)
         call validate_observable_state(s1, current_observables, 'LOCALINIT', 0, temperature_index)

         if (NMESS <= 0) then
            call accumulate_local_nematic_average(s1)
            return
         end if

         do_bind_time = WRITEBINDTIME .and. (iconf <= NTRACECONF) .and. (abs(T - TMIN) < 1.0D-8)
         if (do_bind_time) then
            write (bt_trace_file, '("cfbt",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (77, file=bt_trace_file, status='replace')
            write (77, *) 'iconf T'
            write (77, '(I4,1x,E12.6)') iconf, T
            write (77, *) 'isweep, binderx, bindery'
            disorder_stripe_x_second_moment = 0.D0
            disorder_stripe_x_fourth_moment = 0.D0
            disorder_stripe_y_second_moment = 0.D0
            disorder_stripe_y_fourth_moment = 0.D0
         end if

         do isweep = 1, NMESS
            call metropolis_sweep(s1, beta, flip_counts_this_config(temperature_index, :), current_observables)
            call validate_observable_state(s1, current_observables, 'LOCALMEAS', isweep, temperature_index)
            call accumulate_local_nematic_average(s1)
            call spatial_mean_local_nematic(s1, current_local_nematic_mean)
            sample_magnetization(isweep) = current_observables%magnetization
            sample_energy(isweep) = current_observables%energy
            sample_stripe_x(isweep) = current_observables%stripe_x
            sample_stripe_y(isweep) = current_observables%stripe_y
            sample_local_nematic_mean(isweep) = current_local_nematic_mean

            if (do_bind_time) then
               disorder_stripe_x_second_moment = disorder_stripe_x_second_moment + current_observables%stripe_x**2
               disorder_stripe_x_fourth_moment = disorder_stripe_x_fourth_moment + current_observables%stripe_x**4
               disorder_stripe_y_second_moment = disorder_stripe_y_second_moment + current_observables%stripe_y**2
               disorder_stripe_y_fourth_moment = disorder_stripe_y_fourth_moment + current_observables%stripe_y**4

               mean_stripe_x_second_moment = disorder_stripe_x_second_moment/isweep
               mean_stripe_x_fourth_moment = disorder_stripe_x_fourth_moment/isweep
               mean_stripe_y_second_moment = disorder_stripe_y_second_moment/isweep
               mean_stripe_y_fourth_moment = disorder_stripe_y_fourth_moment/isweep

               if (mean_stripe_x_second_moment > 0.D0) then
                  stripe_x_binder_value = 1.D0 - &
     &               mean_stripe_x_fourth_moment/(3.0D0*mean_stripe_x_second_moment**2)
               else
                  stripe_x_binder_value = 0.D0
               end if
               if (mean_stripe_y_second_moment > 0.D0) then
                  stripe_y_binder_value = 1.D0 - &
     &               mean_stripe_y_fourth_moment/(3.0D0*mean_stripe_y_second_moment**2)
               else
                  stripe_y_binder_value = 0.D0
               end if

               write (77, '(I8,1x,E12.6,1x,E12.6)') isweep, stripe_x_binder_value, stripe_y_binder_value
            end if

         end do
         if (do_bind_time) then
            close (77)
         end if

      end subroutine measurement_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function ising_binder(second_moment, fourth_moment)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b), intent(in) :: second_moment, fourth_moment

         if (second_moment > 0.D0) then
            ising_binder = 1.D0 - fourth_moment/(3.D0*second_moment**2)
         else
            ising_binder = 0.D0
         end if
      end function ising_binder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function two_component_binder(second_moment, fourth_moment)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b), intent(in) :: second_moment, fourth_moment

         if (second_moment > 0.D0) then
            two_component_binder = 2.D0 - fourth_moment/(second_moment**2)
         else
            two_component_binder = 0.D0
         end if
      end function two_component_binder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine store_config_observables(temperature_index, beta_for_temperature, summary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: temperature_index
         real(r8b), intent(in) :: beta_for_temperature
         type(config_observable_summary), intent(in) :: summary
         real(r8b) :: magnetic_susceptibility_value, stripe_susceptibility_value
         real(r8b) :: magnetization_binder_value, stripe_binder_value
         real(r8b) :: stripe_x_binder_value, stripe_y_binder_value
         real(r8b) :: specific_heat_value, local_nematic_binder_value, stripe_difference_value

         ! This is the last step for one disorder config at one temperature.
         ! Everything below feeds the disorder averages and the per-config files.
         stripe_difference_value = summary%mean_abs_stripe_x**2 - summary%mean_abs_stripe_y**2
         magnetic_susceptibility_value = (summary%magnetization_second_moment - summary%mean_abs_magnetization**2)*(L3 - nimp)*beta_for_temperature
         stripe_susceptibility_value = (summary%stripe_second_moment - summary%mean_abs_stripe_order**2)*(L3 - nimp)*beta_for_temperature
         magnetization_binder_value = ising_binder(summary%magnetization_second_moment, summary%magnetization_fourth_moment)
         stripe_binder_value = two_component_binder(summary%stripe_second_moment, summary%stripe_fourth_moment)
         stripe_x_binder_value = ising_binder(summary%stripe_x_second_moment, summary%stripe_x_fourth_moment)
         stripe_y_binder_value = ising_binder(summary%stripe_y_second_moment, summary%stripe_y_fourth_moment)
         local_nematic_binder_value = ising_binder(summary%local_nematic_second_moment, summary%local_nematic_fourth_moment)
         specific_heat_value = (summary%energy_second_moment - summary%mean_energy**2)*(L3 - nimp)*beta_for_temperature**2

         config_abs_magnetization(temperature_index) = config_abs_magnetization(temperature_index) + summary%mean_abs_magnetization
         config_abs_magnetization_squared(temperature_index) = config_abs_magnetization_squared(temperature_index) + summary%mean_abs_magnetization**2
         config_magnetization_second_moment(temperature_index) = config_magnetization_second_moment(temperature_index) + summary%magnetization_second_moment
         config_magnetization_fourth_moment(temperature_index) = config_magnetization_fourth_moment(temperature_index) + summary%magnetization_fourth_moment
         if (summary%mean_abs_magnetization > 0.D0) then
            config_log_abs_magnetization(temperature_index) = config_log_abs_magnetization(temperature_index) + log(summary%mean_abs_magnetization)
         else
            config_log_abs_magnetization(temperature_index) = config_log_abs_magnetization(temperature_index) + log(tiny(1.D0))
         end if
         config_magnetic_susceptibility(temperature_index) = config_magnetic_susceptibility(temperature_index) + magnetic_susceptibility_value
         config_magnetization_binder(temperature_index) = config_magnetization_binder(temperature_index) + magnetization_binder_value
         config_stripe_binder(temperature_index) = config_stripe_binder(temperature_index) + stripe_binder_value
         config_magnetic_susceptibility_squared(temperature_index) = config_magnetic_susceptibility_squared(temperature_index) + magnetic_susceptibility_value**2
         config_magnetization_binder_squared(temperature_index) = config_magnetization_binder_squared(temperature_index) + magnetization_binder_value**2
         config_stripe_binder_squared(temperature_index) = config_stripe_binder_squared(temperature_index) + stripe_binder_value**2
         config_energy(temperature_index) = config_energy(temperature_index) + summary%mean_energy
         config_specific_heat(temperature_index) = config_specific_heat(temperature_index) + specific_heat_value
         config_energy_squared(temperature_index) = config_energy_squared(temperature_index) + summary%mean_energy**2
         config_specific_heat_squared(temperature_index) = config_specific_heat_squared(temperature_index) + specific_heat_value**2
         config_stripe_second_moment(temperature_index) = config_stripe_second_moment(temperature_index) + summary%stripe_second_moment
         config_stripe_fourth_moment(temperature_index) = config_stripe_fourth_moment(temperature_index) + summary%stripe_fourth_moment
         config_abs_stripe_order(temperature_index) = config_abs_stripe_order(temperature_index) + summary%mean_abs_stripe_order
         config_abs_stripe_x(temperature_index) = config_abs_stripe_x(temperature_index) + summary%mean_abs_stripe_x
         ! These are copied once for the current disorder config, then stored in *_by_config arrays later.
         config_stripe_x(temperature_index) = summary%mean_abs_stripe_x
         config_stripe_x_second_moment(temperature_index) = config_stripe_x_second_moment(temperature_index) + summary%stripe_x_second_moment
         config_stripe_x_fourth_moment(temperature_index) = config_stripe_x_fourth_moment(temperature_index) + summary%stripe_x_fourth_moment
         config_stripe_x_binder(temperature_index) = config_stripe_x_binder(temperature_index) + stripe_x_binder_value
         config_stripe_x_binder_squared(temperature_index) = config_stripe_x_binder_squared(temperature_index) + stripe_x_binder_value**2
         config_abs_stripe_y(temperature_index) = config_abs_stripe_y(temperature_index) + summary%mean_abs_stripe_y
         config_stripe_y(temperature_index) = summary%mean_abs_stripe_y
         config_abs_stripe(temperature_index) = summary%mean_abs_stripe_order
         config_stripe_difference(temperature_index) = stripe_difference_value
         config_best_energy(temperature_index) = summary%best_energy
         config_final_energy(temperature_index) = summary%final_energy
         config_stripe_y_second_moment(temperature_index) = config_stripe_y_second_moment(temperature_index) + summary%stripe_y_second_moment
         config_stripe_y_fourth_moment(temperature_index) = config_stripe_y_fourth_moment(temperature_index) + summary%stripe_y_fourth_moment
         config_stripe_y_binder(temperature_index) = config_stripe_y_binder(temperature_index) + stripe_y_binder_value
         config_stripe_y_binder_squared(temperature_index) = config_stripe_y_binder_squared(temperature_index) + stripe_y_binder_value**2
         config_abs_stripe_order_squared(temperature_index) = config_abs_stripe_order_squared(temperature_index) + summary%mean_abs_stripe_order**2
         config_stripe_susceptibility(temperature_index) = config_stripe_susceptibility(temperature_index) + stripe_susceptibility_value
         config_stripe_susceptibility_squared(temperature_index) = config_stripe_susceptibility_squared(temperature_index) + stripe_susceptibility_value**2
         config_nematic_order(temperature_index) = config_nematic_order(temperature_index) + summary%mean_nematic_order
         config_nematic_order_squared(temperature_index) = config_nematic_order_squared(temperature_index) + summary%mean_nematic_order**2
         config_local_nematic_mean(temperature_index) = config_local_nematic_mean(temperature_index) + summary%mean_local_nematic
         config_abs_local_nematic_mean(temperature_index) = config_abs_local_nematic_mean(temperature_index) + summary%mean_abs_local_nematic
         config_local_nematic_mean_squared(temperature_index) = config_local_nematic_mean_squared(temperature_index) + summary%mean_local_nematic**2
         config_local_nematic_second_moment(temperature_index) = config_local_nematic_second_moment(temperature_index) + summary%local_nematic_second_moment
         config_local_nematic_fourth_moment(temperature_index) = config_local_nematic_fourth_moment(temperature_index) + summary%local_nematic_fourth_moment
         config_local_nematic_binder(temperature_index) = config_local_nematic_binder(temperature_index) + local_nematic_binder_value
         config_local_nematic_binder_squared(temperature_index) = config_local_nematic_binder_squared(temperature_index) + local_nematic_binder_value**2
         config_local_nematic_by_temp(temperature_index) = summary%mean_local_nematic
         config_abs_local_nematic_by_temp(temperature_index) = summary%mean_abs_local_nematic
         config_local_nematic_second_by_temp(temperature_index) = summary%local_nematic_second_moment
         config_local_nematic_fourth_by_temp(temperature_index) = summary%local_nematic_fourth_moment
         config_local_nematic_binder_by_temp(temperature_index) = local_nematic_binder_value
      end subroutine store_config_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine compute_config_observables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         real(r8b)       :: measured_magnetization, measured_energy
         real(r8b)       :: mean_abs_magnetization, magnetization_second_moment, magnetization_fourth_moment
         real(r8b)       :: measured_stripe_x, measured_stripe_y, mean_abs_stripe_x, mean_abs_stripe_y
         real(r8b)       :: stripe_x_second_moment, stripe_y_second_moment
         real(r8b)       :: stripe_x_fourth_moment, stripe_y_fourth_moment
         real(r8b)       :: stripe_second_moment, stripe_fourth_moment, mean_abs_stripe_order
         real(r8b)       :: mean_energy, energy_second_moment
         real(r8b)       :: mean_nematic_order
         real(r8b)       :: mean_local_nematic, mean_abs_local_nematic
         real(r8b)       :: local_nematic_second_moment, local_nematic_fourth_moment
         real(r8b)       :: hrange                   ! histogram range padding
         type(observable_state) :: measured_observables
         type(config_observable_summary) :: summary
         real(r8b)       :: measured_nbar
         integer(i4b)    :: energy_histogram_bin_index
         real(r8b)       :: best_energy, final_energy

         mean_abs_magnetization = 0.D0
         magnetization_second_moment = 0.D0
         magnetization_fourth_moment = 0.D0
         mean_energy = 0.D0
         energy_second_moment = 0.D0
         measured_stripe_x = 0.D0
         mean_abs_stripe_x = 0.D0
         stripe_x_second_moment = 0.D0
         stripe_x_fourth_moment = 0.D0
         measured_stripe_y = 0.D0
         mean_abs_stripe_y = 0.D0
         stripe_y_second_moment = 0.D0
         stripe_y_fourth_moment = 0.D0
         stripe_second_moment = 0.D0
         stripe_fourth_moment = 0.D0
         mean_abs_stripe_order = 0.D0
         mean_nematic_order = 0.D0
         mean_local_nematic = 0.D0
         mean_abs_local_nematic = 0.D0
         local_nematic_second_moment = 0.D0
         local_nematic_fourth_moment = 0.D0
         best_energy = 0.D0
         final_energy = 0.D0

         if (NMESS <= 0) then
            ! Some tests use zero measurement sweeps; just measure the current spin state.
            call measure_observables(s1, measured_observables)
            call spatial_mean_local_nematic(s1, measured_nbar)
            call summarize_one_measured_state(measured_observables, measured_nbar, summary)
            if (WRITEHIST .eqv. .true.) then
               config_energy_histogram(temperature_index, :) = 0.D0
            end if
         else
            ! initialize histogram for this temperature and disorder config
            ! set bin ranges once per temperature (first time we see this temperature_index)
            if (WRITEHIST .eqv. .true.) then
               config_energy_histogram(temperature_index, :) = 0.D0
               if (energy_histogram_bin_size(temperature_index) <= 0.D0) then
                  energy_histogram_min(temperature_index) = minval(sample_energy)
                  energy_histogram_max(temperature_index) = maxval(sample_energy)
                  if (energy_histogram_max(temperature_index) > energy_histogram_min(temperature_index)) then
                     hrange = energy_histogram_max(temperature_index) - energy_histogram_min(temperature_index)
                     energy_histogram_min(temperature_index) = energy_histogram_min(temperature_index) - 0.50D0*hrange
                     energy_histogram_max(temperature_index) = energy_histogram_max(temperature_index) + 0.50D0*hrange
                  end if
                  if (energy_histogram_max(temperature_index) <= energy_histogram_min(temperature_index)) then
                     energy_histogram_min(temperature_index) = energy_histogram_min(temperature_index) - 1.0D-6
                     energy_histogram_max(temperature_index) = energy_histogram_max(temperature_index) + 1.0D-6
                  end if
                  energy_histogram_bin_size(temperature_index) = (energy_histogram_max(temperature_index) - energy_histogram_min(temperature_index))/BINS
               end if
            end if

            best_energy = minval(sample_energy)
            final_energy = sample_energy(NMESS)

            ! Walk through the saved measurement sweeps and build the thermal moments.
            do isweep = 1, NMESS
               measured_magnetization = sample_magnetization(isweep)
               measured_stripe_x = sample_stripe_x(isweep)
               mean_abs_stripe_x = mean_abs_stripe_x + abs(measured_stripe_x)
               stripe_x_second_moment = stripe_x_second_moment + measured_stripe_x**2
               stripe_x_fourth_moment = stripe_x_fourth_moment + measured_stripe_x**4
               measured_stripe_y = sample_stripe_y(isweep)
               mean_abs_stripe_y = mean_abs_stripe_y + abs(measured_stripe_y)
               stripe_y_second_moment = stripe_y_second_moment + measured_stripe_y**2
               stripe_y_fourth_moment = stripe_y_fourth_moment + measured_stripe_y**4
               stripe_second_moment = stripe_second_moment + measured_stripe_x**2 + measured_stripe_y**2
               stripe_fourth_moment = stripe_fourth_moment + (measured_stripe_x**2 + measured_stripe_y**2)**2
               mean_abs_stripe_order = mean_abs_stripe_order + sqrt(measured_stripe_x**2 + measured_stripe_y**2)
               mean_abs_magnetization = mean_abs_magnetization + abs(measured_magnetization)
               magnetization_second_moment = magnetization_second_moment + measured_magnetization**2
               magnetization_fourth_moment = magnetization_fourth_moment + measured_magnetization**4
               mean_nematic_order = mean_nematic_order + abs(measured_stripe_y**2 - measured_stripe_x**2)
               measured_nbar = sample_local_nematic_mean(isweep)
               mean_local_nematic = mean_local_nematic + measured_nbar
               mean_abs_local_nematic = mean_abs_local_nematic + abs(measured_nbar)
               local_nematic_second_moment = local_nematic_second_moment + measured_nbar**2
               local_nematic_fourth_moment = local_nematic_fourth_moment + measured_nbar**4

               measured_energy = sample_energy(isweep)
               mean_energy = mean_energy + measured_energy
               energy_second_moment = energy_second_moment + measured_energy**2

               if (WRITEHIST .eqv. .true.) then
                  energy_histogram_bin_index = NINT((measured_energy - energy_histogram_min(temperature_index)) &
     &               /energy_histogram_bin_size(temperature_index))
                  if (energy_histogram_bin_index .ge. 1 .and. energy_histogram_bin_index .le. BINS) then
                     config_energy_histogram(temperature_index, energy_histogram_bin_index) = &
     &                  config_energy_histogram(temperature_index, energy_histogram_bin_index) + 1
                  end if
               end if

            end do

            ! Normalize sums into averages before converting them to observables.
            mean_abs_magnetization = mean_abs_magnetization/NMESS
            magnetization_second_moment = magnetization_second_moment/NMESS
            magnetization_fourth_moment = magnetization_fourth_moment/NMESS
            mean_energy = mean_energy/NMESS
            energy_second_moment = energy_second_moment/NMESS
            stripe_second_moment = stripe_second_moment/NMESS
            stripe_fourth_moment = stripe_fourth_moment/NMESS
            mean_abs_stripe_order = mean_abs_stripe_order/NMESS
            mean_nematic_order = mean_nematic_order/NMESS
            mean_local_nematic = mean_local_nematic/NMESS
            mean_abs_local_nematic = mean_abs_local_nematic/NMESS
            local_nematic_second_moment = local_nematic_second_moment/NMESS
            local_nematic_fourth_moment = local_nematic_fourth_moment/NMESS
            mean_abs_stripe_x = mean_abs_stripe_x/NMESS
            mean_abs_stripe_y = mean_abs_stripe_y/NMESS
            stripe_x_second_moment = stripe_x_second_moment/NMESS
            stripe_x_fourth_moment = stripe_x_fourth_moment/NMESS
            stripe_y_second_moment = stripe_y_second_moment/NMESS
            stripe_y_fourth_moment = stripe_y_fourth_moment/NMESS
         end if

         ! Store this temperature's numbers into the current-config arrays.
         if (NMESS > 0) then
            summary%mean_abs_magnetization = mean_abs_magnetization
            summary%magnetization_second_moment = magnetization_second_moment
            summary%magnetization_fourth_moment = magnetization_fourth_moment
            summary%mean_energy = mean_energy
            summary%energy_second_moment = energy_second_moment
            summary%mean_abs_stripe_x = mean_abs_stripe_x
            summary%mean_abs_stripe_y = mean_abs_stripe_y
            summary%stripe_x_second_moment = stripe_x_second_moment
            summary%stripe_x_fourth_moment = stripe_x_fourth_moment
            summary%stripe_y_second_moment = stripe_y_second_moment
            summary%stripe_y_fourth_moment = stripe_y_fourth_moment
            summary%stripe_second_moment = stripe_second_moment
            summary%stripe_fourth_moment = stripe_fourth_moment
            summary%mean_abs_stripe_order = mean_abs_stripe_order
            summary%mean_nematic_order = mean_nematic_order
            summary%mean_local_nematic = mean_local_nematic
            summary%mean_abs_local_nematic = mean_abs_local_nematic
            summary%local_nematic_second_moment = local_nematic_second_moment
            summary%local_nematic_fourth_moment = local_nematic_fourth_moment
            summary%best_energy = best_energy
            summary%final_energy = final_energy
         end if
         call store_config_observables(temperature_index, beta, summary)

         if (WRITEHIST .eqv. .true. .and. NMESS > 0) then
            config_energy_histogram(temperature_index, :) = config_energy_histogram(temperature_index, :)/NMESS
         end if
         config_flip_counts(temperature_index, :) = flip_counts_this_config(temperature_index, :)
         config_phi(:) = p(:)

      end subroutine compute_config_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine build_sweep_site_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: order_index, swap_index, saved_site_index

         do order_index = 0, L3 - 1
            site_order(order_index) = order_index
         end do

         if (UPDATEORDER == 'RND') then
            ! Generate one fresh Fisher-Yates permutation per sweep from the
            ! existing simulation RNG so randomized sweeps stay seed-reproducible.
            do order_index = L3 - 1, 1, -1
               swap_index = int((order_index + 1)*rkiss05())
               saved_site_index = site_order(order_index)
               site_order(order_index) = site_order(swap_index)
               site_order(swap_index) = saved_site_index
            end do
         end if
      end subroutine build_sweep_site_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine compute_local_field(spins, site, field, local_nematic_field)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Local field for the active Hamiltonian, including the phi coupling.
!!    metropolis_sweep uses this helper so the sequential and randomized sweep
!!    variants share identical acceptance physics.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in)  :: spins(0:L3 - 1)
         integer(i4b), intent(in) :: site
         real(r8b), intent(out)   :: field
         real(r8b), intent(out)   :: local_nematic_field

         field = 0.0D0
         local_nematic_field = 0.0D0
         if (occu(m1(site)) == 1) then
            field = field + spins(m1(site))*J1v*jv(site)
            local_nematic_field = local_nematic_field - spins(m1(site))
            field = field - p(m1(site))*spins(m1(site))
         end if
         if (occu(m2(site)) == 1) then
            field = field + spins(m2(site))*J1v*jv(m2(site))
            local_nematic_field = local_nematic_field - spins(m2(site))
            field = field - p(m2(site))*spins(m2(site))
         end if
         if (occu(m3(site)) == 1) then
            field = field + spins(m3(site))*J1h*jh(site)
            local_nematic_field = local_nematic_field + spins(m3(site))
            field = field + p(m3(site))*spins(m3(site))
         end if
         if (occu(m4(site)) == 1) then
            field = field + spins(m4(site))*J1h*jh(m4(site))
            local_nematic_field = local_nematic_field + spins(m4(site))
            field = field + p(m4(site))*spins(m4(site))
         end if

         if (occu(m5(site)) == 1) field = field + spins(m5(site))*J2
         if (occu(m6(site)) == 1) field = field + spins(m6(site))*J2
         if (occu(m7(site)) == 1) field = field + spins(m7(site))*J2
         if (occu(m8(site)) == 1) field = field + spins(m8(site))*J2

         if (occu(m9(site)) == 1) field = field + spins(m9(site))*J3
         if (occu(m10(site)) == 1) field = field + spins(m10(site))*J3

         field = field + (p(site)*local_nematic_field)
      end subroutine compute_local_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine metropolis_sweep(spins, beta_local, flipcounts, observables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Performs a Metropolis sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(inout) :: spins(0:L3 - 1)
         real(r8b), intent(in)       :: beta_local
         integer(i4b), intent(inout) :: flipcounts(0:L3 - 1)
         type(observable_state), intent(inout) :: observables
         real(r8b)   :: field
         real(r8b)   :: local_nematic_field
         real(r8b)   :: norm, delta_mag, delta_en, delta_strx, delta_stry
         integer(i4b) :: order_index
         integer(i1b)  :: old_spin
         logical(ilog) :: flip_accepted

         call build_sweep_site_order
         norm = real(L3 - nimp, r8b)

         do order_index = 0, L3 - 1
            isite = site_order(order_index)
            if (occu(isite) == 1) then
               ! Evaluate the current site's local energy change.
               call compute_local_field(spins, isite, field, local_nematic_field)
               old_spin = spins(isite)
               delta_mag = -2.D0*real(old_spin, r8b)/norm
               delta_en = 2.D0*real(old_spin, r8b)*field/norm
               delta_strx = delta_mag*stripe_parity_x(isite)
               delta_stry = delta_mag*stripe_parity_y(isite)
               flip_accepted = .false.

               if (old_spin == 1) then
                  if (field .le. 0) then
                     flip_accepted = .true.
                  elseif (rkiss05() .lt. exp(-2.D0*field*beta_local)) then
                     flip_accepted = .true.
                  end if
               else
                  if (field .ge. 0) then
                     flip_accepted = .true.
                  elseif (rkiss05() .lt. exp(2.D0*field*beta_local)) then
                     flip_accepted = .true.
                  end if
               end if

               if (flip_accepted) then
                  spins(isite) = -old_spin
                  flipcounts(isite) = flipcounts(isite) + 1
                  observables%magnetization = observables%magnetization + delta_mag
                  observables%energy = observables%energy + delta_en
                  observables%stripe_x = observables%stripe_x + delta_strx
                  observables%stripe_y = observables%stripe_y + delta_stry
               end if
            end if
         end do

      end subroutine metropolis_sweep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measure_observables(spins, observables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Measures magnetization, energy, and stripe order from the full spin state.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in)  :: spins(0:L3 - 1)
         type(observable_state), intent(out) :: observables
         real(r8b)       :: magnetization_sum, energy_numerator
         real(r8b)       :: stripe_x_sum, stripe_y_sum
         real(r8b)       :: field, local_nematic_field
         integer(i4b)    :: layer, site_index

         energy_numerator = 0.D0
         magnetization_sum = 0.D0
         stripe_x_sum = 0.D0
         stripe_y_sum = 0.D0

         do layer = 0, H - 1
            do isite = 0, L2 - 1
               site_index = isite + L2*layer
               if (occu(site_index) == 1) then
                  magnetization_sum = magnetization_sum + spins(site_index)
                  call compute_local_field(spins, site_index, field, local_nematic_field)
                  energy_numerator = energy_numerator - spins(site_index)*field
                  stripe_x_sum = stripe_x_sum + spins(site_index)*(-1)**mod(isite, L)
                  stripe_y_sum = stripe_y_sum + spins(site_index)*(-1)**(isite/L)
               end if
            end do
         end do

         observables%magnetization = magnetization_sum/(L3 - nimp)
         ! The local field is the same one used in Metropolis acceptance, so
         ! the 1/2 prefactor here reports the full simulated Hamiltonian.
         observables%energy = 0.5D0*energy_numerator/(L3 - nimp)
         observables%stripe_x = stripe_x_sum/(L3 - nimp)
         observables%stripe_y = stripe_y_sum/(L3 - nimp)

      end subroutine measure_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_local_nematic_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Writes the measurement-averaged local nematic order at each site.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)     :: ix, iy, iz                   ! x, y, and z coordinates
         character(len=20) :: filename
         real(r8b)         :: sample_norm, measured_mean_local_nematic, expected_mean_local_nematic
         integer(i4b)      :: occ_count

         filename = 'locnem0000.dat'

         write (filename(7:10), '(I4.4)') temperature_index

         open (8, file=filename, status='replace')
         rewind (8)
         sample_norm = real(max(1, local_nematic_sample_count), r8b)
         occ_count = max(1, L3 - nimp)
         measured_mean_local_nematic = sum(local_nematic_sum_by_site(:))/(sample_norm*real(occ_count, r8b))
         expected_mean_local_nematic = 4.D0*(config_stripe_y(temperature_index)**2 - config_stripe_x(temperature_index)**2)
         write (8, '(A,E24.16)') '# Temperature = ', T
         write (8, '(A,I0)') '# disorder configuration = ', iconf
         write (8, '(A,2(E24.16,1x))') '# stripe magnitudes = ', config_stripe_x(temperature_index), config_stripe_y(temperature_index)
         write (8, '(A,E24.16)') '# expected mean locnem = ', expected_mean_local_nematic
         write (8, '(A,E24.16)') '# measured mean locnem = ', measured_mean_local_nematic
         write (8, '(A)') '# ix iy iz isite occu <s(isite)> <locnem(isite)>'
         do iz = 0, H - 1
            do iy = 0, L - 1
               do ix = 0, L - 1
                  isite = L*iy + ix + L2*iz
                  write (8, '(4(I6,1x),I3,1x,2(E18.10,1x))') ix, iy, iz, isite, occu(isite), &
                     spin_sum_by_site(isite)/sample_norm, local_nematic_sum_by_site(isite)/sample_norm
               end do
            end do
         end do
         close (8)

      end subroutine write_local_nematic_order

      subroutine collect_config_data(active_config)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         logical(ilog), intent(in) :: active_config
         integer(i4b), parameter :: BASEFIELDS = 49
         integer(i4b), parameter :: PTFIELDS = 2*(NREPLICA - 1) + 1
         real(r8b)       :: transdata(BASEFIELDS*NTEMP + (NTEMP*BINS) + (NTEMP*L3) + L3 + 1 + PTFIELDS)
         integer(i4b)    :: cfg_id
         integer(i4b)    :: hist_offset, flips_offset, phi_offset, iconf_offset, pt_offset

         hist_offset = BASEFIELDS*NTEMP
         flips_offset = hist_offset + NTEMP*BINS
         phi_offset = flips_offset + NTEMP*L3
         iconf_offset = phi_offset + L3
         pt_offset = iconf_offset + 1
         transdata(:) = 0.D0

#ifdef PARALLEL
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if (myid .ne. 0) then
            if (active_config) then
               ! Pack this rank's config into one real array so MPI only sends one message.
               transdata(1:NTEMP) = config_abs_magnetization(:)
               transdata(NTEMP + 1:2*NTEMP) = config_magnetization_second_moment(:)
               transdata(2*NTEMP + 1:3*NTEMP) = config_magnetization_fourth_moment(:)
               transdata(3*NTEMP + 1:4*NTEMP) = config_magnetic_susceptibility(:)
               transdata(4*NTEMP + 1:5*NTEMP) = config_magnetic_susceptibility_squared(:)
               transdata(5*NTEMP + 1:6*NTEMP) = config_magnetization_binder(:)
               transdata(6*NTEMP + 1:7*NTEMP) = config_magnetization_binder_squared(:)
               transdata(7*NTEMP + 1:8*NTEMP) = config_abs_magnetization_squared(:)
               transdata(8*NTEMP + 1:9*NTEMP) = config_log_abs_magnetization(:)
               transdata(9*NTEMP + 1:10*NTEMP) = config_energy(:)
               transdata(10*NTEMP + 1:11*NTEMP) = config_energy_squared(:)
               transdata(11*NTEMP + 1:12*NTEMP) = config_specific_heat(:)
               transdata(12*NTEMP + 1:13*NTEMP) = config_specific_heat_squared(:)
               transdata(13*NTEMP + 1:14*NTEMP) = config_abs_stripe_order(:)
               transdata(14*NTEMP + 1:15*NTEMP) = config_abs_stripe_x(:)
               transdata(15*NTEMP + 1:16*NTEMP) = config_abs_stripe_y(:)
               transdata(16*NTEMP + 1:17*NTEMP) = config_stripe_second_moment(:)
               transdata(17*NTEMP + 1:18*NTEMP) = config_stripe_fourth_moment(:)
               transdata(18*NTEMP + 1:19*NTEMP) = config_stripe_susceptibility(:)
               transdata(19*NTEMP + 1:20*NTEMP) = config_stripe_susceptibility_squared(:)
               transdata(20*NTEMP + 1:21*NTEMP) = config_stripe_binder(:)
               transdata(21*NTEMP + 1:22*NTEMP) = config_stripe_binder_squared(:)
               transdata(22*NTEMP + 1:23*NTEMP) = config_nematic_order(:)
               transdata(23*NTEMP + 1:24*NTEMP) = config_nematic_order_squared(:)
               transdata(24*NTEMP + 1:25*NTEMP) = config_abs_stripe_order_squared(:)
               transdata(25*NTEMP + 1:26*NTEMP) = config_stripe_x_second_moment(:)
               transdata(26*NTEMP + 1:27*NTEMP) = config_stripe_x_fourth_moment(:)
               transdata(27*NTEMP + 1:28*NTEMP) = config_stripe_x_binder(:)
               transdata(28*NTEMP + 1:29*NTEMP) = config_stripe_x_binder_squared(:)
               transdata(29*NTEMP + 1:30*NTEMP) = config_stripe_y_second_moment(:)
               transdata(30*NTEMP + 1:31*NTEMP) = config_stripe_y_fourth_moment(:)
               transdata(31*NTEMP + 1:32*NTEMP) = config_stripe_y_binder(:)
               transdata(32*NTEMP + 1:33*NTEMP) = config_stripe_y_binder_squared(:)
               transdata(33*NTEMP + 1:34*NTEMP) = config_abs_stripe(:)
               transdata(34*NTEMP + 1:35*NTEMP) = config_stripe_difference(:)
               transdata(35*NTEMP + 1:36*NTEMP) = config_best_energy(:)
               transdata(36*NTEMP + 1:37*NTEMP) = config_final_energy(:)
               transdata(37*NTEMP + 1:38*NTEMP) = config_local_nematic_mean(:)
               transdata(38*NTEMP + 1:39*NTEMP) = config_local_nematic_mean_squared(:)
               transdata(39*NTEMP + 1:40*NTEMP) = config_local_nematic_second_moment(:)
               transdata(40*NTEMP + 1:41*NTEMP) = config_local_nematic_fourth_moment(:)
               transdata(41*NTEMP + 1:42*NTEMP) = config_local_nematic_binder(:)
               transdata(42*NTEMP + 1:43*NTEMP) = config_local_nematic_binder_squared(:)
               transdata(43*NTEMP + 1:44*NTEMP) = config_local_nematic_by_temp(:)
               transdata(44*NTEMP + 1:45*NTEMP) = config_local_nematic_second_by_temp(:)
               transdata(45*NTEMP + 1:46*NTEMP) = config_local_nematic_fourth_by_temp(:)
               transdata(46*NTEMP + 1:47*NTEMP) = config_local_nematic_binder_by_temp(:)
               transdata(47*NTEMP + 1:48*NTEMP) = config_abs_local_nematic_mean(:)
               transdata(48*NTEMP + 1:49*NTEMP) = config_abs_local_nematic_by_temp(:)

               do temperature_index = 1, NTEMP
                  transdata(hist_offset + 1 + (temperature_index - 1)*BINS:hist_offset + temperature_index*BINS) = config_energy_histogram(temperature_index, :)
               end do

               do temperature_index = 1, NTEMP
                  transdata(flips_offset + 1 + (temperature_index - 1)*L3:flips_offset + temperature_index*L3) = &
                     real(config_flip_counts(temperature_index, :), r8b)
               end do

               transdata(phi_offset + 1:phi_offset + L3) = config_phi(:)
               transdata(iconf_offset + 1) = real(iconf, r8b)
               transdata(pt_offset + 1:pt_offset + (NREPLICA - 1)) = real(config_pt_swap_attempts(:), r8b)
               transdata(pt_offset + NREPLICA:pt_offset + 2*(NREPLICA - 1)) = real(config_pt_swap_accepts(:), r8b)
               transdata(pt_offset + 2*(NREPLICA - 1) + 1) = real(config_pt_visits_both_ends, r8b)
            end if

            call MPI_SEND(transdata, size(transdata), MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
         else
            if (active_config) then
               call add_config_to_disorder_sums(iconf)
            end if

            do source_rank = 1, numprocs - 1
               ! Rank 0 unpacks each worker result and adds it to the same global arrays.
               call MPI_RECV(transdata, size(transdata), MPI_DOUBLE_PRECISION, source_rank, 1, MPI_COMM_WORLD, status, ierr)
               cfg_id = nint(transdata(iconf_offset + 1))
               if (cfg_id < 1 .or. cfg_id > NCONF) cycle

               stripe_by_config(cfg_id, 1, :) = transdata(14*NTEMP + 1:15*NTEMP)
               stripe_by_config(cfg_id, 2, :) = transdata(15*NTEMP + 1:16*NTEMP)
               abs_stripe_by_config(cfg_id, :) = transdata(33*NTEMP + 1:34*NTEMP)
               stripe_difference_by_config(cfg_id, :) = transdata(34*NTEMP + 1:35*NTEMP)
               best_energy_by_config(cfg_id, :) = transdata(35*NTEMP + 1:36*NTEMP)
               final_energy_by_config(cfg_id, :) = transdata(36*NTEMP + 1:37*NTEMP)
               local_nematic_mean_by_config(cfg_id, :) = transdata(43*NTEMP + 1:44*NTEMP)
               abs_local_nematic_mean_by_config(cfg_id, :) = transdata(48*NTEMP + 1:49*NTEMP)
               local_nematic_second_moment_by_config(cfg_id, :) = transdata(44*NTEMP + 1:45*NTEMP)
               local_nematic_fourth_moment_by_config(cfg_id, :) = transdata(45*NTEMP + 1:46*NTEMP)
               local_nematic_binder_by_config(cfg_id, :) = transdata(46*NTEMP + 1:47*NTEMP)

               disorder_abs_magnetization(:) = disorder_abs_magnetization(:) + transdata(1:NTEMP)
               disorder_magnetization_second_moment(:) = disorder_magnetization_second_moment(:) + transdata(NTEMP + 1:2*NTEMP)
               disorder_magnetization_fourth_moment(:) = disorder_magnetization_fourth_moment(:) + transdata(2*NTEMP + 1:3*NTEMP)
               disorder_magnetic_susceptibility(:) = disorder_magnetic_susceptibility(:) + transdata(3*NTEMP + 1:4*NTEMP)
               disorder_magnetic_susceptibility_squared(:) = disorder_magnetic_susceptibility_squared(:) + transdata(4*NTEMP + 1:5*NTEMP)
               disorder_magnetization_binder(:) = disorder_magnetization_binder(:) + transdata(5*NTEMP + 1:6*NTEMP)
               disorder_magnetization_binder_squared(:) = disorder_magnetization_binder_squared(:) + transdata(6*NTEMP + 1:7*NTEMP)
               disorder_abs_magnetization_squared(:) = disorder_abs_magnetization_squared(:) + transdata(7*NTEMP + 1:8*NTEMP)
               disorder_log_abs_magnetization(:) = disorder_log_abs_magnetization(:) + transdata(8*NTEMP + 1:9*NTEMP)
               disorder_energy(:) = disorder_energy(:) + transdata(9*NTEMP + 1:10*NTEMP)
               disorder_energy_squared(:) = disorder_energy_squared(:) + transdata(10*NTEMP + 1:11*NTEMP)
               disorder_specific_heat(:) = disorder_specific_heat(:) + transdata(11*NTEMP + 1:12*NTEMP)
               disorder_specific_heat_squared(:) = disorder_specific_heat_squared(:) + transdata(12*NTEMP + 1:13*NTEMP)
               disorder_abs_stripe_order(:) = disorder_abs_stripe_order(:) + transdata(13*NTEMP + 1:14*NTEMP)
               disorder_abs_stripe_x(:) = disorder_abs_stripe_x(:) + transdata(14*NTEMP + 1:15*NTEMP)
               disorder_abs_stripe_y(:) = disorder_abs_stripe_y(:) + transdata(15*NTEMP + 1:16*NTEMP)
               disorder_stripe_second_moment(:) = disorder_stripe_second_moment(:) + transdata(16*NTEMP + 1:17*NTEMP)
               disorder_stripe_fourth_moment(:) = disorder_stripe_fourth_moment(:) + transdata(17*NTEMP + 1:18*NTEMP)
               disorder_stripe_susceptibility(:) = disorder_stripe_susceptibility(:) + transdata(18*NTEMP + 1:19*NTEMP)
               disorder_stripe_susceptibility_squared(:) = disorder_stripe_susceptibility_squared(:) + transdata(19*NTEMP + 1:20*NTEMP)
               disorder_stripe_binder(:) = disorder_stripe_binder(:) + transdata(20*NTEMP + 1:21*NTEMP)
               disorder_stripe_binder_squared(:) = disorder_stripe_binder_squared(:) + transdata(21*NTEMP + 1:22*NTEMP)
               disorder_nematic_order(:) = disorder_nematic_order(:) + transdata(22*NTEMP + 1:23*NTEMP)
               disorder_nematic_order_squared(:) = disorder_nematic_order_squared(:) + transdata(23*NTEMP + 1:24*NTEMP)
               disorder_abs_stripe_order_squared(:) = disorder_abs_stripe_order_squared(:) + transdata(24*NTEMP + 1:25*NTEMP)
               disorder_stripe_x_second_moment(:) = disorder_stripe_x_second_moment(:) + transdata(25*NTEMP + 1:26*NTEMP)
               disorder_stripe_x_fourth_moment(:) = disorder_stripe_x_fourth_moment(:) + transdata(26*NTEMP + 1:27*NTEMP)
               disorder_stripe_x_binder(:) = disorder_stripe_x_binder(:) + transdata(27*NTEMP + 1:28*NTEMP)
               disorder_stripe_x_binder_squared(:) = disorder_stripe_x_binder_squared(:) + transdata(28*NTEMP + 1:29*NTEMP)
               disorder_stripe_y_second_moment(:) = disorder_stripe_y_second_moment(:) + transdata(29*NTEMP + 1:30*NTEMP)
               disorder_stripe_y_fourth_moment(:) = disorder_stripe_y_fourth_moment(:) + transdata(30*NTEMP + 1:31*NTEMP)
               disorder_stripe_y_binder(:) = disorder_stripe_y_binder(:) + transdata(31*NTEMP + 1:32*NTEMP)
               disorder_stripe_y_binder_squared(:) = disorder_stripe_y_binder_squared(:) + transdata(32*NTEMP + 1:33*NTEMP)
               disorder_local_nematic_mean(:) = disorder_local_nematic_mean(:) + transdata(37*NTEMP + 1:38*NTEMP)
               disorder_abs_local_nematic_mean(:) = disorder_abs_local_nematic_mean(:) + transdata(47*NTEMP + 1:48*NTEMP)
               disorder_local_nematic_mean_squared(:) = disorder_local_nematic_mean_squared(:) + transdata(38*NTEMP + 1:39*NTEMP)
               disorder_local_nematic_second_moment(:) = disorder_local_nematic_second_moment(:) + transdata(39*NTEMP + 1:40*NTEMP)
               disorder_local_nematic_fourth_moment(:) = disorder_local_nematic_fourth_moment(:) + transdata(40*NTEMP + 1:41*NTEMP)
               disorder_local_nematic_binder(:) = disorder_local_nematic_binder(:) + transdata(41*NTEMP + 1:42*NTEMP)
               disorder_local_nematic_binder_squared(:) = disorder_local_nematic_binder_squared(:) + transdata(42*NTEMP + 1:43*NTEMP)
               do temperature_index = 1, NTEMP
                  flip_counts_by_config(cfg_id, temperature_index, :) = int(&
     &               transdata(flips_offset + 1 + (temperature_index - 1)*L3:&
     &               flips_offset + temperature_index*L3)&
     &            )
               end do
               phi_by_config(cfg_id, :) = transdata(phi_offset + 1:phi_offset + L3)

               do temperature_index = 1, NTEMP
                disorder_energy_histogram(temperature_index, :) = disorder_energy_histogram(temperature_index, :) + transdata(hist_offset + (temperature_index - 1)*BINS + 1:hist_offset + temperature_index*BINS)
               end do
               disorder_pt_swap_attempts(:) = disorder_pt_swap_attempts(:) + int(transdata(pt_offset + 1:pt_offset + (NREPLICA - 1)))
               disorder_pt_swap_accepts(:) = disorder_pt_swap_accepts(:) + int(transdata(pt_offset + NREPLICA:pt_offset + 2*(NREPLICA - 1)))
               disorder_pt_visits_both_ends = disorder_pt_visits_both_ends + nint(transdata(pt_offset + 2*(NREPLICA - 1) + 1))

            end do
         end if
#else
         ! Serial mode already has the current config locally, so just add it.
         call add_config_to_disorder_sums(iconf)

#endif
      end subroutine collect_config_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine add_config_to_disorder_sums(config_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: config_id

         ! Save per-config diagnostics first; these become cfxy/cfnd rows.
         stripe_by_config(config_id, 1, :) = config_stripe_x(:)
         stripe_by_config(config_id, 2, :) = config_stripe_y(:)
         abs_stripe_by_config(config_id, :) = config_abs_stripe(:)
         stripe_difference_by_config(config_id, :) = config_stripe_difference(:)
         best_energy_by_config(config_id, :) = config_best_energy(:)
         final_energy_by_config(config_id, :) = config_final_energy(:)
         local_nematic_mean_by_config(config_id, :) = config_local_nematic_by_temp(:)
         abs_local_nematic_mean_by_config(config_id, :) = config_abs_local_nematic_by_temp(:)
         local_nematic_second_moment_by_config(config_id, :) = config_local_nematic_second_by_temp(:)
         local_nematic_fourth_moment_by_config(config_id, :) = config_local_nematic_fourth_by_temp(:)
         local_nematic_binder_by_config(config_id, :) = config_local_nematic_binder_by_temp(:)

         ! Add this config to the disorder-level sums used by the summary files.
         disorder_abs_magnetization(:) = disorder_abs_magnetization(:) + config_abs_magnetization(:)
         disorder_magnetization_second_moment(:) = disorder_magnetization_second_moment(:) + config_magnetization_second_moment(:)
         disorder_magnetization_fourth_moment(:) = disorder_magnetization_fourth_moment(:) + config_magnetization_fourth_moment(:)
         disorder_magnetic_susceptibility(:) = disorder_magnetic_susceptibility(:) + config_magnetic_susceptibility(:)
         disorder_magnetic_susceptibility_squared(:) = disorder_magnetic_susceptibility_squared(:) + config_magnetic_susceptibility_squared(:)
         disorder_magnetization_binder(:) = disorder_magnetization_binder(:) + config_magnetization_binder(:)
         disorder_magnetization_binder_squared(:) = disorder_magnetization_binder_squared(:) + config_magnetization_binder_squared(:)
         disorder_abs_magnetization_squared(:) = disorder_abs_magnetization_squared(:) + config_abs_magnetization_squared(:)
         disorder_log_abs_magnetization(:) = disorder_log_abs_magnetization(:) + config_log_abs_magnetization(:)
         disorder_energy(:) = disorder_energy(:) + config_energy(:)
         disorder_energy_squared(:) = disorder_energy_squared(:) + config_energy_squared(:)
         disorder_specific_heat(:) = disorder_specific_heat(:) + config_specific_heat(:)
         disorder_specific_heat_squared(:) = disorder_specific_heat_squared(:) + config_specific_heat_squared(:)
         disorder_abs_stripe_order(:) = disorder_abs_stripe_order(:) + config_abs_stripe_order(:)
         disorder_abs_stripe_x(:) = disorder_abs_stripe_x(:) + config_abs_stripe_x(:)
         disorder_abs_stripe_y(:) = disorder_abs_stripe_y(:) + config_abs_stripe_y(:)
         disorder_stripe_second_moment(:) = disorder_stripe_second_moment(:) + config_stripe_second_moment(:)
         disorder_stripe_fourth_moment(:) = disorder_stripe_fourth_moment(:) + config_stripe_fourth_moment(:)
         disorder_stripe_susceptibility(:) = disorder_stripe_susceptibility(:) + config_stripe_susceptibility(:)
         disorder_stripe_susceptibility_squared(:) = disorder_stripe_susceptibility_squared(:) + config_stripe_susceptibility_squared(:)
         disorder_stripe_binder(:) = disorder_stripe_binder(:) + config_stripe_binder(:)
         disorder_stripe_binder_squared(:) = disorder_stripe_binder_squared(:) + config_stripe_binder_squared(:)
         disorder_nematic_order(:) = disorder_nematic_order(:) + config_nematic_order(:)
         disorder_nematic_order_squared(:) = disorder_nematic_order_squared(:) + config_nematic_order_squared(:)
         disorder_abs_stripe_order_squared(:) = disorder_abs_stripe_order_squared(:) + config_abs_stripe_order_squared(:)
         disorder_stripe_x_second_moment(:) = disorder_stripe_x_second_moment(:) + config_stripe_x_second_moment(:)
         disorder_stripe_x_fourth_moment(:) = disorder_stripe_x_fourth_moment(:) + config_stripe_x_fourth_moment(:)
         disorder_stripe_x_binder(:) = disorder_stripe_x_binder(:) + config_stripe_x_binder(:)
         disorder_stripe_x_binder_squared(:) = disorder_stripe_x_binder_squared(:) + config_stripe_x_binder_squared(:)
         disorder_stripe_y_second_moment(:) = disorder_stripe_y_second_moment(:) + config_stripe_y_second_moment(:)
         disorder_stripe_y_fourth_moment(:) = disorder_stripe_y_fourth_moment(:) + config_stripe_y_fourth_moment(:)
         disorder_stripe_y_binder(:) = disorder_stripe_y_binder(:) + config_stripe_y_binder(:)
         disorder_stripe_y_binder_squared(:) = disorder_stripe_y_binder_squared(:) + config_stripe_y_binder_squared(:)
         disorder_local_nematic_mean(:) = disorder_local_nematic_mean(:) + config_local_nematic_mean(:)
         disorder_abs_local_nematic_mean(:) = disorder_abs_local_nematic_mean(:) + config_abs_local_nematic_mean(:)
         disorder_local_nematic_mean_squared(:) = disorder_local_nematic_mean_squared(:) + config_local_nematic_mean_squared(:)
         disorder_local_nematic_second_moment(:) = disorder_local_nematic_second_moment(:) + config_local_nematic_second_moment(:)
         disorder_local_nematic_fourth_moment(:) = disorder_local_nematic_fourth_moment(:) + config_local_nematic_fourth_moment(:)
         disorder_local_nematic_binder(:) = disorder_local_nematic_binder(:) + config_local_nematic_binder(:)
         disorder_local_nematic_binder_squared(:) = disorder_local_nematic_binder_squared(:) + config_local_nematic_binder_squared(:)
         do temperature_index = 1, NTEMP
            flip_counts_by_config(config_id, temperature_index, :) = config_flip_counts(temperature_index, :)
         end do
         phi_by_config(config_id, :) = config_phi(:)
         disorder_pt_swap_attempts(:) = disorder_pt_swap_attempts(:) + config_pt_swap_attempts(:)
         disorder_pt_swap_accepts(:) = disorder_pt_swap_accepts(:) + config_pt_swap_accepts(:)
         disorder_pt_visits_both_ends = disorder_pt_visits_both_ends + config_pt_visits_both_ends
         disorder_energy_histogram(:, :) = disorder_energy_histogram(:, :) + config_energy_histogram(:, :)
      end subroutine add_config_to_disorder_sums

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_run_header(unit_number, size_label, include_xi, include_energy_definition, separator)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: unit_number
         character(len=*), intent(in) :: size_label
         logical(ilog), intent(in) :: include_xi, include_energy_definition
         character(len=*), intent(in) :: separator

         write (unit_number, *) 'program J1J2_6          '
         write (unit_number, *) size_label, L
         write (unit_number, *) 'J1                      ', J1
         write (unit_number, *) 'deltaJ1                 ', deltaJ1
         write (unit_number, *) 'J2                      ', J2
         write (unit_number, *) 'J3                      ', J3
         if (include_xi) then
            write (unit_number, *) 'phi, xi                     ', phi, corr_length
         else
            write (unit_number, *) 'phi                     ', phi
         end if
         write (unit_number, *) 'impurity conc.          ', PP
         write (unit_number, *) 'impurity mode           ', IMPMODE
         write (unit_number, *) 'start                   ', START
         write (unit_number, *) 'stripe axis             ', STRAXIS
         write (unit_number, *) 'sampler mode            ', SAMPLERMODE
         write (unit_number, *) 'update order            ', UPDATEORDER
         write (unit_number, *) 'temperature protocol    ', TEMPPROTOCOL
         write (unit_number, *) 'PT start mode           ', PTSTARTMODE
         write (unit_number, *) 'PT ladder mode          ', PTLADDERMODE
         write (unit_number, *) 'PT replicas             ', NREPLICA
         write (unit_number, *) 'PT swap every           ', PTSWAPEVER
         write (unit_number, *) 'PT measure index        ', PTMEASINDEX
         write (unit_number, *) 'disorder configurations ', NCONF
         write (unit_number, *) 'equilibration steps     ', NEQ
         write (unit_number, *) 'measurement steps       ', NMESS
         write (unit_number, *) 'RNG seed                ', IRINIT
         if (PTLADDERMODE == 'SPLIT') then
            write (unit_number, *) 'T ladder midpoint       ', TLADDERMID
            write (unit_number, *) 'T ladder low-T fraction ', TLADDERFRAC
         end if
         if (include_energy_definition) then
            write (unit_number, *) 'energy definition       ', 'exchange + phi'
         end if
         write (unit_number, *) separator
      end subroutine write_run_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_output_files(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count

#ifdef PARALLEL
         if (myid == 0) then
#endif
            ! Each writer owns one output file family.  The formats stay old-style on purpose.
            call write_energy_file(finished_config_count)
            call write_magnetization_file(finished_config_count)
            call write_nematic_summary_file(finished_config_count)
            call write_stripe_summary_file(finished_config_count)
            call write_histogram_file(finished_config_count)
            call write_config_stripe_file(finished_config_count)
            call write_config_nematic_file(finished_config_count)
            call write_flip_count_file(finished_config_count)
            call write_pt_diagnostics_file(finished_config_count)
            call write_phi_table_file(finished_config_count)
#ifdef PARALLEL
         end if
#endif
      end subroutine write_output_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_energy_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count

         open (7, file=cfenfile, status='replace')
         rewind (7)
         call write_run_header(7, 'linear system size      ', .true., .true., '-------------------------')
         write (7, *) 'disorder configurations processed ', finished_config_count
         write (7, *) '   T      [<energy>]      std.en/sqrt(iconf)      [<spec.heat>]      std.spec.heat/sqrt(iconf)   '
         do temperature_index = 1, NTEMP
            T = temperature_for_index(temperature_index)
            write (7, '(1x,25(es24.16,1x))') t, disorder_energy(temperature_index)/finished_config_count, &
               (sqrt(disorder_energy_squared(temperature_index)/finished_config_count - (disorder_energy(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               disorder_specific_heat(temperature_index)/finished_config_count, &
               (sqrt(disorder_specific_heat_squared(temperature_index)/finished_config_count - (disorder_specific_heat(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count)
         end do
         close (7)
      end subroutine write_energy_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_magnetization_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count

         open (7, file=cfmagfile, status='replace')
         rewind (7)
         call write_run_header(7, 'spatial system size     ', .true., .false., '------------------------')
         write (7, *) 'disorder configurations processed ', finished_config_count
         write (7, *) '   T      [<mag>]      [<mag>^2]      std.mag/sqrt(iconf)      [<susc>]', &
            ' .     std.susc/sqrt(iconf)      [<Binder>]      std.bin/sqrt(iconf) ', &
            '[log <mag>]  Global.Binder'
         do temperature_index = 1, NTEMP
            T = temperature_for_index(temperature_index)
            write (7, '(1x,25(e12.6,1x))') t, disorder_abs_magnetization(temperature_index)/finished_config_count, disorder_abs_magnetization_squared(temperature_index)/finished_config_count, &
               (sqrt(disorder_abs_magnetization_squared(temperature_index)/finished_config_count - (disorder_abs_magnetization(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), disorder_magnetic_susceptibility(temperature_index)/finished_config_count, &
               (sqrt(disorder_magnetic_susceptibility_squared(temperature_index)/finished_config_count - (disorder_magnetic_susceptibility(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), disorder_magnetization_binder(temperature_index)/finished_config_count, &
               (sqrt(disorder_magnetization_binder_squared(temperature_index)/finished_config_count - (disorder_magnetization_binder(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), disorder_log_abs_magnetization(temperature_index)/finished_config_count, &
               1 - disorder_magnetization_fourth_moment(temperature_index)/finished_config_count/(3*(disorder_magnetization_second_moment(temperature_index)/finished_config_count)**2)
         end do
         close (7)
      end subroutine write_magnetization_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_nematic_summary_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         real(r8b) :: local_nematic_binder_value, mean_local_nematic_error
         real(r8b) :: mean_abs_config_local_nematic, mean_abs_config_local_nematic_error
         real(r8b) :: abs_local_nematic_error

         open (7, file=cfnmfile, status='replace')
         rewind (7)
         call write_run_header(7, 'spatial system size     ', .true., .false., '------------------------')
         write (7, *) 'disorder configurations processed ', finished_config_count
         write (7, *) '   T      [<n>]      [|<n>|]      [<|n|>]      [<n>].Error      [|<n>|].Error      [<|n|>].Error      [<n^2>]      [<n^4>]      Global.Nematic.Binder      Nematic.Binder      Nematic.Binder.Error'
         do temperature_index = 1, NTEMP
            T = temperature_for_index(temperature_index)
            local_nematic_binder_value = ising_binder(disorder_local_nematic_second_moment(temperature_index)/finished_config_count, disorder_local_nematic_fourth_moment(temperature_index)/finished_config_count)
            mean_local_nematic_error = sqrt(max(0.D0, &
               disorder_local_nematic_mean_squared(temperature_index)/finished_config_count - &
               (disorder_local_nematic_mean(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count)
            mean_abs_config_local_nematic = &
               sum(abs(local_nematic_mean_by_config(1:finished_config_count, temperature_index)))/finished_config_count
            mean_abs_config_local_nematic_error = sqrt(max(0.D0, &
               sum((abs(local_nematic_mean_by_config(1:finished_config_count, temperature_index)) - &
               mean_abs_config_local_nematic)**2)/finished_config_count))/sqrt(1.D0*finished_config_count)
            abs_local_nematic_error = sqrt(max(0.D0, &
               sum((abs_local_nematic_mean_by_config(1:finished_config_count, temperature_index) - &
               disorder_abs_local_nematic_mean(temperature_index)/finished_config_count)**2)/finished_config_count))/ &
               sqrt(1.D0*finished_config_count)
            write (7, '(1x,25(e12.6,1x))') T, &
               disorder_local_nematic_mean(temperature_index)/finished_config_count, &
               mean_abs_config_local_nematic, &
               disorder_abs_local_nematic_mean(temperature_index)/finished_config_count, &
               mean_local_nematic_error, &
               mean_abs_config_local_nematic_error, &
               abs_local_nematic_error, &
               disorder_local_nematic_second_moment(temperature_index)/finished_config_count, &
               disorder_local_nematic_fourth_moment(temperature_index)/finished_config_count, &
               local_nematic_binder_value, &
               disorder_local_nematic_binder(temperature_index)/finished_config_count, &
               (sqrt(disorder_local_nematic_binder_squared(temperature_index)/finished_config_count - (disorder_local_nematic_binder(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count)
         end do
         close (7)
      end subroutine write_nematic_summary_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_stripe_summary_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         integer(i4b) :: config_id
         real(r8b) :: coexprod_sum, coexprod2_sum, coexangle_sum, coexangle2_sum

         open (7, file=cfstrfile, status='replace')
         rewind (7)
         call write_run_header(7, 'spatial system size     ', .true., .false., '------------------------')
         write (7, *) 'disorder configurations processed ', finished_config_count
         write (7, *) '   T      absstr    str2      str4      nem      Global.Strip.Binder      Strip.susc ', &
            ' Local.Strip.Binder     Local.Strip.Binder.Error      Strip.OP.error      stri.susc.error', &
            ' strx   strx2   strx4   Local.Strx.Binder    Local.Strx.Binder.Error ', &
            ' stry   stry2   stry4   Local.Stry.Binder    Local.Stry.Binder.Error', &
            ' absstrx_absstry  absstrx_absstry.Error  coex_angle  coex_angle.Error'
         do temperature_index = 1, NTEMP
            T = temperature_for_index(temperature_index)
            coexprod_sum = 0.D0
            coexprod2_sum = 0.D0
            coexangle_sum = 0.D0
            coexangle2_sum = 0.D0
            do config_id = 1, finished_config_count
               coexprod_sum = coexprod_sum + stripe_by_config(config_id, 1, temperature_index)*stripe_by_config(config_id, 2, temperature_index)
               coexangle_sum = coexangle_sum + stripe_coexistence_angle(stripe_by_config(config_id, 1, temperature_index), stripe_by_config(config_id, 2, temperature_index))
               coexprod2_sum = coexprod2_sum + (stripe_by_config(config_id, 1, temperature_index)*stripe_by_config(config_id, 2, temperature_index))**2
               coexangle2_sum = coexangle2_sum + stripe_coexistence_angle(stripe_by_config(config_id, 1, temperature_index), stripe_by_config(config_id, 2, temperature_index))**2
            end do
            write (7, '(1x,25(e12.6,1x))') T, &
               disorder_abs_stripe_order(temperature_index)/finished_config_count, &
               disorder_stripe_second_moment(temperature_index)/finished_config_count, &
               disorder_stripe_fourth_moment(temperature_index)/finished_config_count, &
               disorder_nematic_order(temperature_index)/finished_config_count, &
               two_component_binder(disorder_stripe_second_moment(temperature_index)/finished_config_count, disorder_stripe_fourth_moment(temperature_index)/finished_config_count), &
               disorder_stripe_susceptibility(temperature_index)/finished_config_count, &
               disorder_stripe_binder(temperature_index)/finished_config_count, &
               (sqrt(disorder_stripe_binder_squared(temperature_index)/finished_config_count - (disorder_stripe_binder(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               (sqrt(disorder_abs_stripe_order_squared(temperature_index)/finished_config_count - (disorder_abs_stripe_order(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               (sqrt(disorder_stripe_susceptibility_squared(temperature_index)/finished_config_count - (disorder_stripe_susceptibility(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               disorder_abs_stripe_x(temperature_index)/finished_config_count, &
               disorder_stripe_x_second_moment(temperature_index)/finished_config_count, &
               disorder_stripe_x_fourth_moment(temperature_index)/finished_config_count, &
               disorder_stripe_x_binder(temperature_index)/finished_config_count, &
               (sqrt(disorder_stripe_x_binder_squared(temperature_index)/finished_config_count - (disorder_stripe_x_binder(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               disorder_abs_stripe_y(temperature_index)/finished_config_count, &
               disorder_stripe_y_second_moment(temperature_index)/finished_config_count, &
               disorder_stripe_y_fourth_moment(temperature_index)/finished_config_count, &
               disorder_stripe_y_binder(temperature_index)/finished_config_count, &
               (sqrt(disorder_stripe_y_binder_squared(temperature_index)/finished_config_count - (disorder_stripe_y_binder(temperature_index)/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               coexprod_sum/finished_config_count, &
               (sqrt(coexprod2_sum/finished_config_count - (coexprod_sum/finished_config_count)**2))/sqrt(1.D0*finished_config_count), &
               coexangle_sum/finished_config_count, &
               (sqrt(coexangle2_sum/finished_config_count - (coexangle_sum/finished_config_count)**2))/sqrt(1.D0*finished_config_count)
         end do
         close (7)
      end subroutine write_stripe_summary_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_histogram_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count

         if (.not. WRITEHIST) return
         open (7, file=cfhsfile, status='replace')
         rewind (7)
         call write_run_header(7, 'linear system size      ', .false., .true., '-------------------------')
         write (7, *) 'disorder configurations processed ', finished_config_count
         write (7, *) 'T  BIN   [<energy>]  freq'
         do temperature_index = 1, NTEMP
            T = temperature_for_index(temperature_index)
            do ibin = 1, BINS - 1
               write (7, '(1x, 1(e12.6, 1x), I4, 1x, 2(e12.6, 1x))') T, ibin, &
                  energy_histogram_min(temperature_index) + energy_histogram_bin_size(temperature_index)*(ibin - 1), &
                  disorder_energy_histogram(temperature_index, ibin)/(finished_config_count)
            end do
         end do
         close (7)
      end subroutine write_histogram_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_config_stripe_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         integer(i4b) :: config_id

         if (.not. WRITESTRIP) return
         open (7, file=cfxyfile, status='replace')
         rewind (7)
         write (7, *) 'iconf, T, strx, stry, absstr, D, best_en, final_en, absstrx_absstry, coex_angle'
         do config_id = 1, finished_config_count
            do temperature_index = 1, NTEMP
               T = temperature_for_index(temperature_index)
               write (7, '(I4,1x,9(E12.6,1x))') config_id, T, stripe_by_config(config_id, 1, temperature_index), stripe_by_config(config_id, 2, temperature_index), &
                  abs_stripe_by_config(config_id, temperature_index), stripe_difference_by_config(config_id, temperature_index), best_energy_by_config(config_id, temperature_index), final_energy_by_config(config_id, temperature_index), &
                  stripe_by_config(config_id, 1, temperature_index)*stripe_by_config(config_id, 2, temperature_index), &
                  stripe_coexistence_angle(stripe_by_config(config_id, 1, temperature_index), stripe_by_config(config_id, 2, temperature_index))
            end do
         end do
         close (7)
      end subroutine write_config_stripe_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_config_nematic_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         integer(i4b) :: config_id

         open (7, file=cfndfile, status='replace')
         rewind (7)
         write (7, *) 'iconf, T, nbar, abs_nbar, nbar2, nbar4, nbinder'
         do config_id = 1, finished_config_count
            do temperature_index = 1, NTEMP
               T = temperature_for_index(temperature_index)
               write (7, '(I4,1x,6(E12.6,1x))') config_id, T, local_nematic_mean_by_config(config_id, temperature_index), &
                  abs_local_nematic_mean_by_config(config_id, temperature_index), local_nematic_second_moment_by_config(config_id, temperature_index), local_nematic_fourth_moment_by_config(config_id, temperature_index), local_nematic_binder_by_config(config_id, temperature_index)
            end do
         end do
         close (7)
      end subroutine write_config_nematic_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_flip_count_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         integer(i4b) :: config_id

         if (.not. WRITENUMFLIPS) return
         open (7, file=cfnumfile, status='replace')
         rewind (7)
         write (7, *) 'iconf, T, isite, numflips'
         do config_id = 1, finished_config_count
            do temperature_index = 1, NTEMP
               T = temperature_for_index(temperature_index)
               do isite = 0, L3 - 1
                  write (7, '(I4,1x,E12.6,1x,I8,1x,I8)') config_id, T, isite, flip_counts_by_config(config_id, temperature_index, isite)
               end do
            end do
         end do
         close (7)
      end subroutine write_flip_count_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_pt_diagnostics_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count

         if (SAMPLERMODE /= 'PT') return
         open (7, file=cfptfile, status='replace')
         rewind (7)
         write (7, *) 'pair, T_low, T_high, attempts, accepts, acceptance'
         do temperature_index = 1, NREPLICA - 1
            if (disorder_pt_swap_attempts(temperature_index) > 0) then
               write (7, '(I4,1x,2(E12.6,1x),2(I12,1x),E12.6)') temperature_index, temp_ladder(temperature_index), temp_ladder(temperature_index + 1), &
                  disorder_pt_swap_attempts(temperature_index), disorder_pt_swap_accepts(temperature_index), &
                  real(disorder_pt_swap_accepts(temperature_index), r8b)/real(disorder_pt_swap_attempts(temperature_index), r8b)
            else
               write (7, '(I4,1x,2(E12.6,1x),2(I12,1x),E12.6)') temperature_index, temp_ladder(temperature_index), temp_ladder(temperature_index + 1), &
                  0, 0, 0.D0
            end if
         end do
         write (7, *) 'replicas_visiting_both_extremes ', disorder_pt_visits_both_ends, &
            real(disorder_pt_visits_both_ends, r8b)/real(max(1, finished_config_count*NREPLICA), r8b)
         close (7)
      end subroutine write_pt_diagnostics_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_phi_table_file(finished_config_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finished_config_count
         integer(i4b) :: config_id

         if (.not. WRITEPHITABLE) return
         open (7, file=cfptable, status='replace')
         rewind (7)
         write (7, *) 'iconf, isite, p'
         do config_id = 1, finished_config_count
            do isite = 0, L3 - 1
               write (7, '(I4,1x,I8,1x,E12.6)') config_id, isite, phi_by_config(config_id, isite)
            end do
         end do
         close (7)
      end subroutine write_phi_table_file

      end program J1J2_6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123
!
!
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
!
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
!
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
!        v0.93    Aug 13, 2012    changed integer representation test to avoid data statements
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION rkiss05()
         implicit none
         integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)   ! 8-byte reals
         integer, parameter      :: i4b = SELECTED_INT_KIND(8)            ! 4-byte integers

         real(r8b), parameter    :: am = 4.656612873077392578d-10       ! multiplier 1/2^31

         real(r8b)             :: rkiss05
         integer(i4b)          :: kiss
         integer(i4b)          :: x, y, z, w              ! working variables for the four generators
         common/kisscom/x, y, z, w

         x = 69069*x + 1327217885
         y = ieor(y, ishft(y, 13)); y = ieor(y, ishft(y, -17)); y = ieor(y, ishft(y, 5))
         z = 18000*iand(z, 65535) + ishft(z, -16)
         w = 30903*iand(w, 65535) + ishft(w, -16)
         kiss = ishft(x + y + ishft(z, 16) + w, -1)
         rkiss05 = kiss*am
      END FUNCTION rkiss05

      SUBROUTINE kissinit(iinit)
         implicit none
         integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)   ! 8-byte reals
         integer, parameter      :: i4b = SELECTED_INT_KIND(8)            ! 4-byte integers

         integer(i4b) idum, ia, im, iq, ir, iinit
         integer(i4b) k, x, y, z, w, c1
         real(r8b) rkiss05, rdum
         parameter(ia=16807, im=2147483647, iq=127773, ir=2836)
         common/kisscom/x, y, z, w

      !!! Test integer representation !!!
         c1 = -8
         c1 = ishftc(c1, -3)
         if (c1 .ne. 536870911) then
            print *, 'Nonstandard integer representation. Stopped.'
            stop
         end if

         idum = iinit
         idum = abs(1099087573*idum)               ! 32-bit LCG to shuffle seeds
         if (idum .eq. 0) idum = 1
         if (idum .ge. IM) idum = IM - 1

         k = (idum)/IQ
         idum = IA*(idum - k*IQ) - IR*k
         if (idum .lt. 0) idum = idum + IM
         if (idum .lt. 1) then
            x = idum + 1
         else
            x = idum
         end if
         k = (idum)/IQ
         idum = IA*(idum - k*IQ) - IR*k
         if (idum .lt. 0) idum = idum + IM
         if (idum .lt. 1) then
            y = idum + 1
         else
            y = idum
         end if
         k = (idum)/IQ
         idum = IA*(idum - k*IQ) - IR*k
         if (idum .lt. 0) idum = idum + IM
         if (idum .lt. 1) then
            z = idum + 1
         else
            z = idum
         end if
         k = (idum)/IQ
         idum = IA*(idum - k*IQ) - IR*k
         if (idum .lt. 0) idum = idum + IM
         if (idum .lt. 1) then
            w = idum + 1
         else
            w = idum
         end if

         rdum = rkiss05()

         return
      end subroutine kissinit
