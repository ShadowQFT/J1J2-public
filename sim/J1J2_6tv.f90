program J1J2_6
   USE, INTRINSIC :: IEEE_EXCEPTIONS
   USE, INTRINSIC :: IEEE_ARITHMETIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Monte Carlo simulation of diluted 3d J1J2 model
!     using Wolff cluster and Metropolis algorithms
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer(i4b), parameter    :: L = 10                               ! linear system size in space and time
   integer(i4b), parameter    :: H = 10                               ! number of layers
   real(r8b), parameter    :: J1 = 1.0D0                           ! nn interaction
   real(r8b), parameter    :: deltaJ1 = 0.0D0                      ! difference of nn interaction
   real(r8b), parameter    :: J1h = J1 + deltaJ1                   ! harizontal J of nn interaction
   real(r8b), parameter    :: J1v = J1 - deltaJ1                   ! vertical J of nn interaction
   real(r8b), parameter    :: J2 = -1.0D0                  ! nnn interaction in terms of nn interaction
   real(r8b), parameter    :: J3 = 1.0D0                           ! interaction between layers
   real(r8b), parameter    :: phi = 2.0D0                         ! domain strength
   real(r8b), parameter    :: corr_length = 3.0D0                  ! correlation length of random field
   real(r8b), parameter    :: TMAX = 14.0D0                          ! temperature loop control - maximum
   real(r8b), parameter    :: TMIN = 8.0D0                          ! temperature loop control - minimum
   real(r8b), parameter    :: DT = 2.0D0                          ! temperature step for LOCAL mode; PT derives spacing from TMIN/TMAX and NREPLICA
   real(r8b), parameter    :: PP = 0.0D0                          ! impurity concentration
   real(r8b), parameter    :: PPh = 0.0D0                            ! harizontal J impurity concentration
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
   logical(ilog), parameter    :: VALIDATE_INCREMENTAL = .false.         ! cross-check maintained observables against full Measurement

   integer(i4b), parameter    :: BINS = 100                             ! number of bins for histogram
   integer(i4b), parameter    :: NTRACECONF = 5                         ! number of disorder configs for low-T binder trace output

   integer(i4b), parameter    :: IRINIT = 1                             ! Random number seed, must be positive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants - do not touch !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer(i4b), parameter    :: L2 = L*L                                    ! 2D system size
   integer(i4b), parameter    :: L3 = L2*H                                   ! 3D system size
   integer(i4b), parameter    :: NREPLICA = 1                                ! PT uses one replica per ladder temperature
   integer(i4b), parameter    :: NTEMP = merge(NREPLICA, 1 + NINT((TMAX - TMIN)/DT), SAMPLERMODE == 'PT')
   integer(i4b), parameter    :: PTMEASINDEX = 1                             ! target temperature index for PT diagnostics (1 = TMIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variable declarations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real(r8b)      :: strxylist(NCONF, 2, NTEMP)                        ! list of stripex and stripey for each disorder config
   real(r8b)      :: cfg_strx(NTEMP), cfg_stry(NTEMP)                  ! current-config stripe averages by temperature
   real(r8b)      :: cfg_absstr_diag(NTEMP), cfg_stripediff(NTEMP)     ! per-config stripe diagnostics by temperature
   real(r8b)      :: cfg_besten(NTEMP), cfg_finalen(NTEMP)             ! per-config best/final energies by temperature
   real(r8b)      :: absstrlist(NCONF, NTEMP), stripedifflist(NCONF, NTEMP)
   real(r8b)      :: bestenlist(NCONF, NTEMP), finalenlist(NCONF, NTEMP)

   integer(i4b)    :: confnumflips(NTEMP, 0:L3 - 1)   ! number of flips per site at each temperature
   integer(i4b)    :: totalnumflips(NCONF, NTEMP, 0:L3 - 1)    ! number of flips per site in each disorder config
   integer(i4b)    :: cfg_numflips(NTEMP, 0:L3 - 1)             ! current-config flip counts by temperature

   real(r8b)       :: confhist(NTEMP, BINS)                           ! list of energy histogram values
   real(r8b)       :: globalhist(NTEMP, BINS)
   real(r8b)       :: hist_en_min_t(NTEMP), hist_en_max_t(NTEMP), hist_bin_size_t(NTEMP)

   integer(i1b)    :: s1(0:L3 - 1)                                      ! Ising variable
   integer(i1b)    :: s1_rep(0:L3 - 1, NREPLICA)                         ! PT replica spin states
   integer(i1b)    :: jh(0:L3 - 1)                                      ! horizontal J validity
   integer(i1b)    :: jv(0:L3 - 1)                                      ! vertical J validity
   integer(i1b)    :: occu(0:L3 - 1)                                    ! occu=1 site occupied, occu=0 site empty
   real(r8b)       :: maglist(NMESS), enlist(NMESS)                    ! list of magnetization and energy for each sweep
   real(r8b)       :: strxlist(NMESS), strylist(NMESS)                 ! list of stripex and stripey for each sweep
   real(r8b)       :: nbarlist(NMESS)                                  ! list of spatially averaged local nematic per sweep

   real(r8b)       :: confmag(NTEMP), confmag2(NTEMP), confmag4(NTEMP)  ! configuration averages of magnetization
   real(r8b)       :: confabsstr(NTEMP), conf2str(NTEMP)
   real(r8b)       :: confstr2(NTEMP), confstr4(NTEMP)                 ! configuration averages of strip
   real(r8b)       :: confstrx(NTEMP), confstry(NTEMP)                ! configuration averages of stripx, stripy
   real(r8b)       :: confstrx2(NTEMP), confstry2(NTEMP)             ! configuration averages of stripx, stripy squared
   real(r8b)       :: confstrx4(NTEMP), confstry4(NTEMP)             ! configuration averages of stripx, stripy to the fourth power
   real(r8b)       :: confsusc(NTEMP), conf2susc(NTEMP)               ! configuration averages of susceptibility
   real(r8b)       :: confsuscstr(NTEMP), conf2suscstr(NTEMP)         ! configuration averages of strip susceptibility
   real(r8b)       :: confUL(NTEMP), conf2UL(NTEMP)                   ! configuration averages of Binder parameter
   real(r8b)       :: confUS(NTEMP), conf2US(NTEMP)                   ! configuration averages of Binder parameter
   real(r8b)       :: confUSX(NTEMP), conf2USX(NTEMP)                 ! configuration averages of stripx binder
   real(r8b)       :: confUSY(NTEMP), conf2USY(NTEMP)                 ! configuration averages of stripy binder
   real(r8b)       :: confen(NTEMP), conf2en(NTEMP)                   ! configuration averages of energy
   real(r8b)       :: confsph(NTEMP), conf2sph(NTEMP)                  ! configuration averages of specific heat
   real(r8b)       :: conf2mag(NTEMP), conflogmag(NTEMP)               !
   real(r8b)       :: confnem(NTEMP), conf2nem(NTEMP)                 ! configuration averages of nematic OP
   real(r8b)       :: confnbar(NTEMP), conf2nbar(NTEMP)               ! configuration averages of spatial mean local nematic
   real(r8b)       :: confnbar2(NTEMP), confnbar4(NTEMP)              ! thermal moments of spatial mean local nematic
   real(r8b)       :: confUBN(NTEMP), conf2UBN(NTEMP)                 ! Binder from spatial mean local nematic

   real(r8b)       :: summag(NTEMP), summag2(NTEMP), summag4(NTEMP)     ! configuration averages of magnetization
   real(r8b)       :: sumabsstr(NTEMP), sum2str(NTEMP)
   real(r8b)       :: sumstr2(NTEMP), sumstr4(NTEMP)                   ! configuration averages of strip
   real(r8b)       :: sumstrx(NTEMP), sumstry(NTEMP)                  ! configuration averages of stripx, stripy
   real(r8b)       :: sumstrx2(NTEMP), sumstry2(NTEMP)             ! configuration averages of stripx, stripy squared
   real(r8b)       :: sumstrx4(NTEMP), sumstry4(NTEMP)             ! configuration averages of stripx, stripy to the fourth power
   real(r8b)       :: sumsusc(NTEMP), sum2susc(NTEMP)                 ! configuration averages of susceptibility
   real(r8b)       :: sumsuscstr(NTEMP), sum2suscstr(NTEMP)           ! configuration averages of strip susceptibility
   real(r8b)       :: sumUL(NTEMP), sum2UL(NTEMP)                     ! configuration averages of Binder parameter
   real(r8b)       :: sumUS(NTEMP), sum2US(NTEMP)                     ! configuration averages of Binder parameter
   real(r8b)       :: sumUSX(NTEMP), sum2USX(NTEMP)                 ! configuration averages of stripx binder
   real(r8b)       :: sumUSY(NTEMP), sum2USY(NTEMP)                 ! configuration averages of stripy binder
   real(r8b)       :: sumen(NTEMP), sum2en(NTEMP)                     ! configuration averages of energy
   real(r8b)       :: sumsph(NTEMP), sum2sph(NTEMP)                    ! configuration averages of specific heat
   real(r8b)       :: sum2mag(NTEMP), sumlogmag(NTEMP)
   real(r8b)       :: sumnem(NTEMP), sum2nem(NTEMP)
   real(r8b)       :: sumnbar(NTEMP), sum2nbar(NTEMP)
   real(r8b)       :: sumnbar2(NTEMP), sumnbar4(NTEMP)
   real(r8b)       :: sumUBN(NTEMP), sum2UBN(NTEMP)

   real(r8b)       :: nbarcfglist(NCONF, NTEMP), nbar2cfglist(NCONF, NTEMP)
   real(r8b)       :: nbar4cfglist(NCONF, NTEMP), nbindcfglist(NCONF, NTEMP)

   real(r8b)       :: avclsize              ! average cluster size in Wolff algorithm
   real(r8b)       :: padd                  ! addition probability for Wolff cluster

   real(r8b)       :: p(0:L3 - 1)             ! random field at each site
   real(r8b)       :: jointp(NCONF, 0:L3 - 1)        ! joined list of random fields over all disorder configs
   real(r8b)       :: cfg_phi(0:L3 - 1)               ! current-config random field values

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
   integer(i4b)    :: isite                 !  site, time indices        !!! isite = L*iy+ix
   integer(i4b)    :: layer                 !  layer index

   real(r8b)       :: T, beta               ! Monte Carlo temperature
   real(r8b)       :: temp_ladder(NTEMP), beta_ladder(NTEMP)            ! explicit ladder used by PT and output helpers
   real(r8b)       :: locspin_accum(0:L3 - 1), locnem_accum(0:L3 - 1)   ! measurement-averaged site diagnostics
   real(r8b)       :: stripe_parity_x(0:L3 - 1), stripe_parity_y(0:L3 - 1)
   real(r8b)       :: cfg_nbar(NTEMP), cfg_nbar2(NTEMP), cfg_nbar4(NTEMP), cfg_nbind(NTEMP)
   integer(i4b)    :: locnem_samples
   integer(i4b)    :: itemp
   integer(i4b)    :: site_order(0:L3 - 1)
   integer(i4b)    :: rep_at_temp(NREPLICA), temp_of_rep(NREPLICA)      ! PT ownership maps between configs and ladder slots
   integer(i4b)    :: cfg_ptswap_attempts(NREPLICA - 1), cfg_ptswap_accepts(NREPLICA - 1)
   integer(i4b)    :: sum_ptswap_attempts(NREPLICA - 1), sum_ptswap_accepts(NREPLICA - 1)
   integer(i4b)    :: cfg_ptvisitboth, sum_ptvisitboth

   integer(i4b)    :: init
   integer(kind=8) :: config_clock_start, collect_clock_start
   logical(ilog)   :: active_config

   real(r8b), external     :: rkiss05
   external kissinit

! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
   include 'mpif.h'
   integer(i4b)              :: ierr
   integer(i4b)              :: id, myid                  ! CPU index
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

!      print *,'J1J2_6'

#ifdef PARALLEL
   call Setup_mpi                                  ! Start MPI
   call Setup_pt_parallel
#endif

   call Setup_neighbor_table                       ! build neighbor table for square lattice with helical BC
   call Setup_stripe_parities
   call Setup_conf_averages                        ! zero the arrays containing configuration averages
   call Validate_run_options
   call Setup_temperature_ladder

#ifdef PARALLEL
   if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
      max_cfg_per_rank = (NCONF + pt_num_groups - 1)/pt_num_groups
   else
      max_cfg_per_rank = (NCONF + numprocs - 1)/numprocs
   end if
   disorder_loop: do confstep = 0, max_cfg_per_rank - 1
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
      disorder_loop: do iconf = 1, NCONF
         active_config = .true.
         print *, 'dis. conf.', iconf
#endif
         if (.not. active_config) then
#ifdef PARALLEL
            call Collect_data(active_config)
#endif
            cycle disorder_loop
         end if
         call Reset_current_config_averages
         call Timer_start(config_clock_start)
         call Write_progress('CONFIG_START', 0, 0)
         init = IRINIT + iconf - 1                             ! initialize RNG
         call kissinit(init)
         call Setup_vacancies                            ! Determine random positions of vacancies
         call Mark_vacancy_flipcounts                  ! Mark vacancy sites so they stay out of flip-count diagnostics
         if (SAMPLERMODE == 'LOCAL' .and. TEMPPROTOCOL == 'CARRY') then
            call Initialize_spins                        ! Carried protocol preserves the original single initialization per disorder config.
         end if
         call Setup_phi_table                            ! Set the random phi values for each site

         ! swap above 2 lines to ensure same initial config for different impurity configs

! Loop over temperatures
         if (SAMPLERMODE == 'PT') then
            ! PT handles the whole ladder for this disorder sample in one pass.
            if (PTPARALLELMODE == 'GROUPED') then
#ifdef PARALLEL
               call Grouped_pt_run
               if (WRITELOCNEM .eqv. .true. .and. iconf == 1 .and. pt_group_master) call Local_nematic_OP
#else
               print *, 'ERROR: PTPARALLELMODE=GROUPED requires PARALLEL build.'
               stop 1
#endif
            else
               call Initialize_pt_state
               call Pt_equilibration_loop
               call Pt_measurement_loop
               if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call Snapshot_pt_target_replica
               if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call Local_nematic_OP
            end if
         else
            temperature_loop: do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               if (myid == 0) then
                  print *, 'Running T= ', T, 'itemp= ', itemp
               end if
               call Prepare_temperature_state
               beta = 1.D0/T
               confnumflips(itemp, :) = 0
               where (occu(:) == 0)                ! mark vacancies so they don't contribute to flip-count diagnostics
                  confnumflips(itemp, :) = -1
               end where

               call Equilibration_loop
               call Measurement_loop
               call Process_data

               if (WRITELOCNEM .eqv. .true. .and. iconf == 1) call Local_nematic_OP

            end do temperature_loop
         end if

         

         !if(WRITELOCNEM.eqv..true.) call Domain_energy

         !print *, 'collecting data for disorder config ', iconf
         call Timer_start(collect_clock_start)
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
#ifdef PARALLEL
            call Collect_data(active_config .and. pt_group_master)
#else
            call Collect_data(active_config)
#endif
         else
            call Collect_data(active_config)
         end if
         call Write_timing('COLLECT_DATA', collect_clock_start)
         call Write_timing('CONFIG_TOTAL', config_clock_start)
         !print *, "num flips per site:", totalnumflips(:,:)

#ifdef PARALLEL
         if (SAMPLERMODE == 'PT' .and. PTPARALLELMODE == 'GROUPED') then
            completed_config_count = min(NCONF, (confstep + 1)*pt_num_groups)
         else
            completed_config_count = min(NCONF, (confstep + 1)*numprocs)
         end if
#else
         completed_config_count = iconf
#endif
         call Write_files(completed_config_count)

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

#ifdef PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

      end subroutine Setup_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_pt_parallel
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

      end subroutine Setup_pt_parallel
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Timer_start(timer_count)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(kind=8), intent(out) :: timer_count

         call system_clock(timer_count)
      end subroutine Timer_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Write_timing(label, timer_count)
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
         call flush(6)
      end subroutine Write_timing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Write_progress(phase, current, total)
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
         close(98)

         write (*, '(A,1x,A,1x,A,I0,1x,A,I0,1x,A,I0)') 'progress', trim(phase), 'iconf', iconf, &
     &      'step', current, 'of', total
         call flush(6)
      end subroutine Write_progress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_neighbor_table
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

      end subroutine Setup_neighbor_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_phi_table
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

         open (99, file='phi_sites.dat', status='replace')
         write (99, *) '# ix iy layer isite phi'

         do layer = 0, H - 1
            do iy = 0, L - 1
               do ix = 0, L - 1
                  isite = L*iy + ix + L2*layer
                  p(isite) = field(ix, iy, layer)
                  write (99, '(3I8,1X,E20.12)') ix, iy, layer, field(ix, iy, layer)
                  !print *, 'phi_site=', p(isite)
               end do
            end do
         end do

         close (99)

      end subroutine Setup_phi_table


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_conf_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         strxylist(:, :, :) = 0.D0
         cfg_strx(:) = 0.D0
         cfg_stry(:) = 0.D0
         cfg_absstr_diag(:) = 0.D0
         cfg_stripediff(:) = 0.D0
         cfg_besten(:) = 0.D0
         cfg_finalen(:) = 0.D0
         cfg_nbar(:) = 0.D0
         cfg_nbar2(:) = 0.D0
         cfg_nbar4(:) = 0.D0
         cfg_nbind(:) = 0.D0
         absstrlist(:, :) = 0.D0
         stripedifflist(:, :) = 0.D0
         bestenlist(:, :) = 0.D0
         finalenlist(:, :) = 0.D0
         nbarcfglist(:, :) = 0.D0
         nbar2cfglist(:, :) = 0.D0
         nbar4cfglist(:, :) = 0.D0
         nbindcfglist(:, :) = 0.D0

         confnumflips(:, :) = 0
         totalnumflips(:, :, :) = 0
         cfg_numflips(:, :) = 0

         confmag(:) = 0.D0
         confmag2(:) = 0.D0
         confmag4(:) = 0.D0
         confsusc(:) = 0.D0
         confUL(:) = 0.D0
         conf2susc(:) = 0.D0
         conf2UL(:) = 0.D0
         confen(:) = 0.D0
         confsph(:) = 0.D0
         conf2en(:) = 0.D0
         conf2sph(:) = 0.D0
         confUS(:) = 0.D0
         conf2US(:) = 0.D0
         confstr2(:) = 0.D0
         confstr4(:) = 0.D0
         confabsstr(:) = 0.D0
         confstrx(:) = 0.D0
         confstry(:) = 0.D0
         confstrx2(:) = 0.D0
         confstry2(:) = 0.D0
         confstrx4(:) = 0.D0
         confstry4(:) = 0.D0
         conf2str(:) = 0.D0
         confsuscstr(:) = 0.D0
         conf2suscstr(:) = 0.D0
         conflogmag(:) = 0.D0
         confnem(:) = 0.D0
         conf2nem(:) = 0.D0
         confUSX(:) = 0.D0
         conf2USX(:) = 0.D0
         confUSY(:) = 0.D0
         conf2USY(:) = 0.D0
         confnbar(:) = 0.D0
         conf2nbar(:) = 0.D0
         confnbar2(:) = 0.D0
         confnbar4(:) = 0.D0
         confUBN(:) = 0.D0
         conf2UBN(:) = 0.D0

         summag(:) = 0.D0
         summag2(:) = 0.D0
         summag4(:) = 0.D0
         sumsusc(:) = 0.D0
         sum2susc(:) = 0.D0
         sumUL(:) = 0.D0
         sum2UL(:) = 0.D0
         sum2mag(:) = 0.D0
         sumlogmag(:) = 0.D0
         sumen(:) = 0.D0
         sum2en(:) = 0.D0
         sumsph(:) = 0.D0
         sum2sph(:) = 0.D0
         sumabsstr(:) = 0.D0
         sumstrx(:) = 0.D0
         sumstry(:) = 0.D0
         sumstr2(:) = 0.D0
         sumstr4(:) = 0.D0
         sumsuscstr(:) = 0.D0
         sum2suscstr(:) = 0.D0
         sumUS(:) = 0.D0
         sum2US(:) = 0.D0
         sumnem(:) = 0.D0
         sum2nem(:) = 0.D0
         sum2str(:) = 0.D0
         sumstrx2(:) = 0.D0
         sumstrx4(:) = 0.D0
         sumUSX(:) = 0.D0
         sum2USX(:) = 0.D0
         sumstry2(:) = 0.D0
         sumstry4(:) = 0.D0
         sumUSY(:) = 0.D0
         sum2USY(:) = 0.D0
         sumnbar(:) = 0.D0
         sum2nbar(:) = 0.D0
         sumnbar2(:) = 0.D0
         sumnbar4(:) = 0.D0
         sumUBN(:) = 0.D0
         sum2UBN(:) = 0.D0

         confhist(:, :) = 0.D0
         globalhist(:, :) = 0.D0
         hist_en_min_t(:) = 0.D0
         hist_en_max_t(:) = 0.D0
         hist_bin_size_t(:) = 0.D0
         cfg_phi(:) = 0.D0
         s1_rep(:, :) = 1
         rep_at_temp(:) = 0
         temp_of_rep(:) = 0
         cfg_ptswap_attempts(:) = 0
         cfg_ptswap_accepts(:) = 0
         sum_ptswap_attempts(:) = 0
         sum_ptswap_accepts(:) = 0
         cfg_ptvisitboth = 0
         sum_ptvisitboth = 0
         temp_ladder(:) = 0.D0
         beta_ladder(:) = 0.D0
         locspin_accum(:) = 0.D0
         locnem_accum(:) = 0.D0
         locnem_samples = 0

      end subroutine Setup_conf_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Reset_current_config_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         cfg_strx(:) = 0.D0
         cfg_stry(:) = 0.D0
         cfg_absstr_diag(:) = 0.D0
         cfg_stripediff(:) = 0.D0
         cfg_besten(:) = 0.D0
         cfg_finalen(:) = 0.D0
         cfg_nbar(:) = 0.D0
         cfg_nbar2(:) = 0.D0
         cfg_nbar4(:) = 0.D0
         cfg_nbind(:) = 0.D0

         confnumflips(:, :) = 0
         cfg_numflips(:, :) = 0

         confmag(:) = 0.D0
         confmag2(:) = 0.D0
         confmag4(:) = 0.D0
         confsusc(:) = 0.D0
         confUL(:) = 0.D0
         conf2susc(:) = 0.D0
         conf2UL(:) = 0.D0
         confen(:) = 0.D0
         confsph(:) = 0.D0
         conf2en(:) = 0.D0
         conf2sph(:) = 0.D0
         confUS(:) = 0.D0
         conf2US(:) = 0.D0
         confstr2(:) = 0.D0
         confstr4(:) = 0.D0
         confabsstr(:) = 0.D0
         confstrx(:) = 0.D0
         confstry(:) = 0.D0
         confstrx2(:) = 0.D0
         confstry2(:) = 0.D0
         confstrx4(:) = 0.D0
         confstry4(:) = 0.D0
         conf2str(:) = 0.D0
         confsuscstr(:) = 0.D0
         conf2suscstr(:) = 0.D0
         conflogmag(:) = 0.D0
         confnem(:) = 0.D0
         conf2nem(:) = 0.D0
         confUSX(:) = 0.D0
         conf2USX(:) = 0.D0
         confUSY(:) = 0.D0
         conf2USY(:) = 0.D0
         confnbar(:) = 0.D0
         conf2nbar(:) = 0.D0
         confnbar2(:) = 0.D0
         confnbar4(:) = 0.D0
         confUBN(:) = 0.D0
         conf2UBN(:) = 0.D0

         confhist(:, :) = 0.D0
         cfg_phi(:) = 0.D0
         cfg_ptswap_attempts(:) = 0
         cfg_ptswap_accepts(:) = 0
         cfg_ptvisitboth = 0
         locspin_accum(:) = 0.D0
         locnem_accum(:) = 0.D0
         locnem_samples = 0

      end subroutine Reset_current_config_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Validate_run_options
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
      end subroutine Validate_run_options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_temperature_ladder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: it, ninterval, nlow, nhigh

         ninterval = NTEMP - 1
         if (SAMPLERMODE == 'PT' .and. PTLADDERMODE == 'SPLIT' .and. ninterval > 0) then
            nlow = nint(ninterval*TLADDERFRAC)
            nlow = max(1, min(ninterval - 1, nlow))
            nhigh = ninterval - nlow

            temp_ladder(1) = TMIN
            do it = 1, nlow
               temp_ladder(it + 1) = TMIN + (TLADDERMID - TMIN)*real(it, r8b)/real(nlow, r8b)
            end do
            do it = 1, nhigh
               temp_ladder(nlow + 1 + it) = TLADDERMID + (TMAX - TLADDERMID)*real(it, r8b)/real(nhigh, r8b)
            end do
         elseif (SAMPLERMODE == 'PT') then
            if (NTEMP <= 1) then
               temp_ladder(1) = TMIN
            else
               do it = 1, NTEMP
                  temp_ladder(it) = TMIN + (TMAX - TMIN)*real(it - 1, r8b)/real(NTEMP - 1, r8b)
               end do
            end if
         else
            do it = 1, NTEMP
               temp_ladder(it) = TMIN + (it - 1)*DT
            end do
         end if

         do it = 1, NTEMP
            beta_ladder(it) = 1.D0/temp_ladder(it)
         end do
      end subroutine Setup_temperature_ladder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_stripe_parities
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
      end subroutine Setup_stripe_parities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function Temperature_for_index(it)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: it

         if (SAMPLERMODE == 'PT') then
            Temperature_for_index = temp_ladder(it)
         elseif (START == 'HOT') then
            Temperature_for_index = temp_ladder(NTEMP - it + 1)
         else
            Temperature_for_index = temp_ladder(it)
         end if
      end function Temperature_for_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Reset_local_nematic_average
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         locspin_accum(:) = 0.D0
         locnem_accum(:) = 0.D0
         locnem_samples = 0
      end subroutine Reset_local_nematic_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Accumulate_local_nematic_average(spins)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         real(r8b)                :: locnem_site

         do isite = 0, L3 - 1
            if (occu(isite) == 1) then
               locspin_accum(isite) = locspin_accum(isite) + real(spins(isite), r8b)
               locnem_site = real(spins(isite)*spins(m3(isite))*occu(m3(isite)), r8b) - &
                  real(spins(isite)*spins(m1(isite))*occu(m1(isite)), r8b) + &
                  real(spins(isite)*spins(m4(isite))*occu(m4(isite)), r8b) - &
                  real(spins(isite)*spins(m2(isite))*occu(m2(isite)), r8b)
               locnem_accum(isite) = locnem_accum(isite) + locnem_site
            end if
         end do
         locnem_samples = locnem_samples + 1
      end subroutine Accumulate_local_nematic_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Spatial_mean_local_nematic(spins, mean_locnem)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         real(r8b), intent(out)   :: mean_locnem
         integer(i4b)             :: occ_count
         real(r8b)                :: locnem_site, locnem_sum

         locnem_sum = 0.D0
         occ_count = max(1, L3 - nimp)

         do isite = 0, L3 - 1
            if (occu(isite) == 0) cycle
            locnem_site = real(spins(isite)*spins(m3(isite))*occu(m3(isite)), r8b) - &
               real(spins(isite)*spins(m1(isite))*occu(m1(isite)), r8b) + &
               real(spins(isite)*spins(m4(isite))*occu(m4(isite)), r8b) - &
               real(spins(isite)*spins(m2(isite))*occu(m2(isite)), r8b)
            locnem_sum = locnem_sum + locnem_site
         end do

         mean_locnem = locnem_sum/real(occ_count, r8b)
      end subroutine Spatial_mean_local_nematic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b) function Stripe_coexistence_angle(strx, stry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b), intent(in) :: strx, stry
         real(r8b)             :: denom

         denom = strx**2 + stry**2
         if (denom > 0.D0) then
            ! This equals 0.5*(1-cos(4 theta)): 1 on the diagonal and 0 on either axis.
            Stripe_coexistence_angle = 4.D0*strx**2*stry**2/(denom**2)
         else
            Stripe_coexistence_angle = 0.D0
         end if
      end function Stripe_coexistence_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Initialize_observables(spins, curmag, curen, curstrx, curstry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         real(r8b), intent(out)   :: curmag, curen, curstrx, curstry

         call Measurement(spins, curmag, curen, curstrx, curstry)
      end subroutine Initialize_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Validate_observable_state(spins, curmag, curen, curstrx, curstry, context, idx1, idx2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in) :: spins(0:L3 - 1)
         real(r8b), intent(in)    :: curmag, curen, curstrx, curstry
         character(len=*), intent(in) :: context
         integer(i4b), intent(in) :: idx1, idx2
         real(r8b)                :: refmag, refen, refstrx, refstry
         real(r8b), parameter     :: tol = 1.0D-10

         if (.not. VALIDATE_INCREMENTAL) return

         call Measurement(spins, refmag, refen, refstrx, refstry)
         if (abs(curmag - refmag) > tol .or. abs(curen - refen) > tol .or. &
     &       abs(curstrx - refstrx) > tol .or. abs(curstry - refstry) > tol) then
            print *, 'ERROR: incremental observable mismatch in ', trim(context), idx1, idx2
            print *, ' current=', curmag, curen, curstrx, curstry
            print *, ' reference=', refmag, refen, refstrx, refstry
            stop 1
         end if
      end subroutine Validate_observable_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Prepare_temperature_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: temp_seed

         if (TEMPPROTOCOL == 'INDEP') then
            ! Reset the RNG before each temperature so fresh-start tests are reproducible
            ! and independent of how many random numbers earlier temperatures consumed.
            temp_seed = IRINIT + (iconf - 1)*NTEMP + itemp
            call kissinit(temp_seed)
            call Initialize_spins
         end if
      end subroutine Prepare_temperature_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Setup_vacancies
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
                  print *, 'adding impurity'
                  nimp = nimp + 1
                  print *, nimp
               end if
               print *, 'one cycle passed'
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

      end subroutine Setup_vacancies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Mark_vacancy_flipcounts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: temp_slot

         do temp_slot = 1, NTEMP
            where (occu(:) == 0)
               confnumflips(temp_slot, :) = -1
            end where
         end do
      end subroutine Mark_vacancy_flipcounts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Initialize_spin_array(spins)
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

      end subroutine Initialize_spin_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Initialize_spins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call Initialize_spin_array(s1)
      end subroutine Initialize_spins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Initialize_pt_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: rep
         integer(i4b) :: replica_seed

         ! Reset PT diagnostics for the new disorder sample.
         cfg_ptswap_attempts(:) = 0
         cfg_ptswap_accepts(:) = 0
         cfg_ptvisitboth = 0
         confhist(:, :) = 0.D0
         hist_en_min_t(:) = 0.D0
         hist_en_max_t(:) = 0.D0
         hist_bin_size_t(:) = 0.D0
         confnumflips(:, :) = 0
         call Mark_vacancy_flipcounts

         ! Start with replica k sitting at ladder slot k.
         do rep = 1, NREPLICA
            rep_at_temp(rep) = rep
            temp_of_rep(rep) = rep
            replica_seed = IRINIT + (iconf - 1)*NREPLICA + rep
            call kissinit(replica_seed)
            ! Each replica gets its own spin configuration on the same disorder sample.
            call Initialize_spin_array(s1_rep(0:L3 - 1, rep))
         end do

         ! Use a dedicated post-initialization seed so PT dynamics remain
         ! reproducible even when the number of replicas changes later.
         call kissinit(IRINIT + 100000 + iconf)
      end subroutine Initialize_pt_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pt_equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep
         integer(i4b) :: swap_phase
         integer(i4b) :: rep
         integer(i4b) :: progress_every
         integer(kind=8) :: phase_clock_start
         real(r8b)    :: rep_mag(NREPLICA), rep_en(NREPLICA), rep_strx(NREPLICA), rep_stry(NREPLICA)

         call Timer_start(phase_clock_start)
         do rep = 1, NREPLICA
            call Initialize_observables(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep))
            call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep), 'PTEQINIT', 0, rep)
         end do

         ! Alternate odd and even neighbor pairs on successive swap passes.
         swap_phase = 1
         progress_every = max(1, NEQ/10)
         call Write_progress('PT_EQUIL', 0, NEQ)
         do isweep = 1, NEQ
            call Pt_cycle(swap_phase, rep_mag, rep_en, rep_strx, rep_stry)
            if (swap_phase == 1) then
               swap_phase = 2
            else
               swap_phase = 1
            end if
            if (mod(isweep, progress_every) == 0 .or. isweep == NEQ) call Write_progress('PT_EQUIL', isweep, NEQ)
         end do
         call Write_timing('PT_EQUIL', phase_clock_start)
      end subroutine Pt_equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pt_cycle(swap_phase, rep_mag, rep_en, rep_strx, rep_stry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: swap_phase
         real(r8b), intent(inout) :: rep_mag(NREPLICA), rep_en(NREPLICA), rep_strx(NREPLICA), rep_stry(NREPLICA)
         integer(i4b)             :: slot, rep, isweep

         ! Evolve each ladder slot at its current temperature before swapping.
         do slot = 1, NREPLICA
            rep = rep_at_temp(slot)
            do isweep = 1, PTSWAPEVER
               call Metro_sweep(s1_rep(0:L3 - 1, rep), beta_ladder(slot), confnumflips(slot, :), &
     &            rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep))
            end do
            call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep), 'PTCYCLE', slot, rep)
         end do

         ! Then try the active set of nearest-neighbor replica exchanges.
         call Pt_attempt_swaps(swap_phase, rep_en)
      end subroutine Pt_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pt_attempt_swaps(swap_phase, rep_en)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: swap_phase
         real(r8b), intent(in)    :: rep_en(NREPLICA)
         integer(i4b)             :: slot, slot_start, rep_lo, rep_hi, tmp
         real(r8b)                :: delta
         real(r8b)                :: total_en_lo, total_en_hi

         if (swap_phase == 1) then
            slot_start = 1
         else
            slot_start = 2
         end if

         ! Swap neighboring ladder slots only, alternating odd/even bonds.
         do slot = slot_start, NREPLICA - 1, 2
            rep_lo = rep_at_temp(slot)
            rep_hi = rep_at_temp(slot + 1)
            cfg_ptswap_attempts(slot) = cfg_ptswap_attempts(slot) + 1

            total_en_lo = rep_en(rep_lo)*(L3 - nimp)
            total_en_hi = rep_en(rep_hi)*(L3 - nimp)
            delta = (beta_ladder(slot) - beta_ladder(slot + 1))*(total_en_lo - total_en_hi)

            ! Accept with the standard replica-exchange Metropolis rule.
            if (delta >= 0.D0 .or. rkiss05() < exp(delta)) then
               cfg_ptswap_accepts(slot) = cfg_ptswap_accepts(slot) + 1
               tmp = rep_at_temp(slot)
               rep_at_temp(slot) = rep_at_temp(slot + 1)
               rep_at_temp(slot + 1) = tmp
               ! Keep the inverse map in sync after the swap.
               temp_of_rep(rep_lo) = slot + 1
               temp_of_rep(rep_hi) = slot
            end if
         end do
      end subroutine Pt_attempt_swaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pt_measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, slot, rep, nmeasure, histindex
         integer(i4b) :: swap_phase
         integer(i4b) :: progress_every
         integer(kind=8) :: phase_clock_start
         real(r8b)    :: rep_mag(NREPLICA), rep_en(NREPLICA), rep_strx(NREPLICA), rep_stry(NREPLICA), rep_nbar(NREPLICA)
         real(r8b)    :: mag_sum(NTEMP), absmag_sum(NTEMP), mag2_sum(NTEMP), mag4_sum(NTEMP)
         real(r8b)    :: en_sum(NTEMP), en2_sum(NTEMP), strx_abs_sum(NTEMP), stry_abs_sum(NTEMP)
         real(r8b)    :: strx2_sum(NTEMP), strx4_sum(NTEMP), stry2_sum(NTEMP), stry4_sum(NTEMP)
         real(r8b)    :: str2_sum(NTEMP), str4_sum(NTEMP), absstr_sum(NTEMP), nem_sum(NTEMP)
         real(r8b)    :: nbar_sum(NTEMP), nbar2_sum(NTEMP), nbar4_sum(NTEMP)
         real(r8b)    :: best_en(NTEMP), final_en(NTEMP), stripediff(NTEMP)
         real(r8b)    :: amag, aen, astripx, astripy, anbar, str2, absstr, nem
         real(r8b)    :: mag, absmag, mag2, mag4, en, en2, totalstrx, totalstry
         real(r8b)    :: strx2, strx4, stry2, stry4, str4, suscabs, susstrip, UL, US, USX, USY, sph
         real(r8b)    :: nbar, nbar2, nbar4, UBN
         real(r8b)    :: target_sumstrx2, target_sumstrx4, target_sumstry2, target_sumstry4
         real(r8b)    :: avgstrx2, avgstrx4, avgstry2, avgstry4, bindx, bindy
         logical(ilog) :: do_bind_time, visited_cold(NREPLICA), visited_hot(NREPLICA), do_trace
         character(len=21) :: bt_trace_file, trace_file

         mag_sum(:) = 0.D0
         absmag_sum(:) = 0.D0
         mag2_sum(:) = 0.D0
         mag4_sum(:) = 0.D0
         en_sum(:) = 0.D0
         en2_sum(:) = 0.D0
         strx_abs_sum(:) = 0.D0
         stry_abs_sum(:) = 0.D0
         strx2_sum(:) = 0.D0
         strx4_sum(:) = 0.D0
         stry2_sum(:) = 0.D0
         stry4_sum(:) = 0.D0
         str2_sum(:) = 0.D0
         str4_sum(:) = 0.D0
         absstr_sum(:) = 0.D0
         nem_sum(:) = 0.D0
         nbar_sum(:) = 0.D0
         nbar2_sum(:) = 0.D0
         nbar4_sum(:) = 0.D0
         stripediff(:) = 0.D0
         best_en(:) = huge(1.D0)
         final_en(:) = 0.D0
         cfg_phi(:) = p(:)
         cfg_numflips(:, :) = confnumflips(:, :)
         call Reset_local_nematic_average
         call Timer_start(phase_clock_start)

         ! Track whether a replica has visited both ends of the ladder.
         visited_cold(:) = .false.
         visited_hot(:) = .false.
         do_bind_time = WRITEBINDTIME .and. (iconf <= NTRACECONF)
         do_trace = (iconf <= NTRACECONF)

         if (do_bind_time) then
            write (bt_trace_file, '("cfbt",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (77, file=bt_trace_file, status='replace')
            write (77, *) 'iconf T'
            write (77, '(I4,1x,E12.6)') iconf, temp_ladder(PTMEASINDEX)
            write (77, *) 'isweep, binderx, bindery'
            target_sumstrx2 = 0.D0
            target_sumstrx4 = 0.D0
            target_sumstry2 = 0.D0
            target_sumstry4 = 0.D0
         end if

         if (do_trace) then
            write (trace_file, '("cftr",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (78, file=trace_file, status='replace')
            write (78, *) 'isweep replica temp_index T'
         end if

         if (NMESS <= 0) then
            call Write_progress('PT_MEAS', 0, 0)
            do rep = 1, NREPLICA
               call Initialize_observables(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep))
               call Spatial_mean_local_nematic(s1_rep(0:L3 - 1, rep), rep_nbar(rep))
            end do
            do rep = 1, NREPLICA
               call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep), 'PTINIT', 0, rep)
            end do
            call Accumulate_local_nematic_average(s1_rep(0:L3 - 1, rep_at_temp(PTMEASINDEX)))
            nmeasure = 1
         else
            ! In PT mode NEQ and NMESS count whole PT cycles, not raw sweeps.
            nmeasure = NMESS
            swap_phase = 1
            progress_every = max(1, NMESS/10)
            do rep = 1, NREPLICA
               call Initialize_observables(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep))
               call Spatial_mean_local_nematic(s1_rep(0:L3 - 1, rep), rep_nbar(rep))
               call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag(rep), rep_en(rep), rep_strx(rep), rep_stry(rep), 'PTINIT', 0, rep)
            end do
            call Write_progress('PT_MEAS', 0, NMESS)
            do isweep = 1, NMESS
               call Pt_cycle(swap_phase, rep_mag, rep_en, rep_strx, rep_stry)

               do rep = 1, NREPLICA
                  if (temp_of_rep(rep) == 1) visited_cold(rep) = .true.
                  if (temp_of_rep(rep) == NREPLICA) visited_hot(rep) = .true.
                  call Spatial_mean_local_nematic(s1_rep(0:L3 - 1, rep), rep_nbar(rep))
               end do

               if (do_trace) then
                  do rep = 1, min(NREPLICA, 3)
                     write (78, '(I8,1x,I4,1x,I4,1x,E12.6)') isweep, rep, temp_of_rep(rep), temp_ladder(temp_of_rep(rep))
                  end do
               end if

               if (do_bind_time) then
                  rep = rep_at_temp(PTMEASINDEX)
                  target_sumstrx2 = target_sumstrx2 + rep_strx(rep)**2
                  target_sumstrx4 = target_sumstrx4 + rep_strx(rep)**4
                  target_sumstry2 = target_sumstry2 + rep_stry(rep)**2
                  target_sumstry4 = target_sumstry4 + rep_stry(rep)**4

                  avgstrx2 = target_sumstrx2/isweep
                  avgstrx4 = target_sumstrx4/isweep
                  avgstry2 = target_sumstry2/isweep
                  avgstry4 = target_sumstry4/isweep

                  if (avgstrx2 > 0.D0) then
                     bindx = 1.D0 - avgstrx4/(3.D0*avgstrx2**2)
                  else
                     bindx = 0.D0
                  end if
                  if (avgstry2 > 0.D0) then
                     bindy = 1.D0 - avgstry4/(3.D0*avgstry2**2)
                  else
                     bindy = 0.D0
                  end if
                  write (77, '(I8,1x,E12.6,1x,E12.6)') isweep, bindx, bindy
               end if

               call Accumulate_local_nematic_average(s1_rep(0:L3 - 1, rep_at_temp(PTMEASINDEX)))

               ! Accumulate observables by temperature slot after swaps.
               do slot = 1, NTEMP
                  rep = rep_at_temp(slot)
                  amag = rep_mag(rep)
                  aen = rep_en(rep)
                  astripx = rep_strx(rep)
                  astripy = rep_stry(rep)
                  anbar = rep_nbar(rep)
                  str2 = astripx**2 + astripy**2
                  absstr = sqrt(str2)
                  nem = abs(astripy**2 - astripx**2)

                  mag_sum(slot) = mag_sum(slot) + amag
                  absmag_sum(slot) = absmag_sum(slot) + abs(amag)
                  mag2_sum(slot) = mag2_sum(slot) + amag**2
                  mag4_sum(slot) = mag4_sum(slot) + amag**4
                  en_sum(slot) = en_sum(slot) + aen
                  en2_sum(slot) = en2_sum(slot) + aen**2
                  strx_abs_sum(slot) = strx_abs_sum(slot) + abs(astripx)
                  stry_abs_sum(slot) = stry_abs_sum(slot) + abs(astripy)
                  strx2_sum(slot) = strx2_sum(slot) + astripx**2
                  strx4_sum(slot) = strx4_sum(slot) + astripx**4
                  stry2_sum(slot) = stry2_sum(slot) + astripy**2
                  stry4_sum(slot) = stry4_sum(slot) + astripy**4
                  str2_sum(slot) = str2_sum(slot) + str2
                  str4_sum(slot) = str4_sum(slot) + str2**2
                  absstr_sum(slot) = absstr_sum(slot) + absstr
                  nem_sum(slot) = nem_sum(slot) + nem
                  nbar_sum(slot) = nbar_sum(slot) + anbar
                  nbar2_sum(slot) = nbar2_sum(slot) + anbar**2
                  nbar4_sum(slot) = nbar4_sum(slot) + anbar**4
                  best_en(slot) = min(best_en(slot), aen)
                  final_en(slot) = aen

                  if (WRITEHIST .eqv. .true.) then
                     if (hist_bin_size_t(slot) <= 0.D0) then
                        hist_en_min_t(slot) = aen - (abs(aen) + 1.D0)
                        hist_en_max_t(slot) = aen + (abs(aen) + 1.D0)
                        hist_bin_size_t(slot) = (hist_en_max_t(slot) - hist_en_min_t(slot))/BINS
                     end if
                     histindex = 1 + int((aen - hist_en_min_t(slot))/hist_bin_size_t(slot))
                     histindex = max(1, min(BINS, histindex))
                     confhist(slot, histindex) = confhist(slot, histindex) + 1.D0
                  end if
               end do

               if (swap_phase == 1) then
                  swap_phase = 2
               else
                  swap_phase = 1
               end if
               if (mod(isweep, progress_every) == 0 .or. isweep == NMESS) call Write_progress('PT_MEAS', isweep, NMESS)
            end do
         end if
         call Write_timing('PT_MEAS', phase_clock_start)

         if (do_bind_time) close(77)
         if (do_trace) close(78)

         ! Count replicas that completed at least one cold-to-hot visitation.
         cfg_ptvisitboth = count(visited_cold(:) .and. visited_hot(:))

         do slot = 1, NTEMP
            rep = rep_at_temp(slot)
            if (NMESS <= 0) then
               amag = rep_mag(rep)
               aen = rep_en(rep)
               astripx = rep_strx(rep)
               astripy = rep_stry(rep)
               mag = abs(amag)
               absmag = abs(amag)
               mag2 = amag**2
               mag4 = amag**4
               en = aen
               en2 = aen**2
               totalstrx = abs(astripx)
               totalstry = abs(astripy)
               strx2 = astripx**2
               strx4 = astripx**4
               stry2 = astripy**2
               stry4 = astripy**4
               str2 = strx2 + stry2
               str4 = str2**2
               absstr = sqrt(str2)
               nem = abs(astripy**2 - astripx**2)
               nbar = rep_nbar(rep)
               nbar2 = nbar**2
               nbar4 = nbar**4
               best_en(slot) = aen
               final_en(slot) = aen
            else
               mag = abs(mag_sum(slot))/nmeasure
               absmag = absmag_sum(slot)/nmeasure
               mag2 = mag2_sum(slot)/nmeasure
               mag4 = mag4_sum(slot)/nmeasure
               en = en_sum(slot)/nmeasure
               en2 = en2_sum(slot)/nmeasure
               totalstrx = strx_abs_sum(slot)/nmeasure
               totalstry = stry_abs_sum(slot)/nmeasure
               strx2 = strx2_sum(slot)/nmeasure
               strx4 = strx4_sum(slot)/nmeasure
               stry2 = stry2_sum(slot)/nmeasure
               stry4 = stry4_sum(slot)/nmeasure
               str2 = str2_sum(slot)/nmeasure
               str4 = str4_sum(slot)/nmeasure
               absstr = absstr_sum(slot)/nmeasure
               nem = nem_sum(slot)/nmeasure
               nbar = nbar_sum(slot)/nmeasure
               nbar2 = nbar2_sum(slot)/nmeasure
               nbar4 = nbar4_sum(slot)/nmeasure
            end if

            stripediff(slot) = totalstrx**2 - totalstry**2
            suscabs = (mag2 - absmag**2)*(L3 - nimp)*beta_ladder(slot)
            susstrip = (str2 - absstr**2)*(L3 - nimp)*beta_ladder(slot)
            if (mag2 > 0.D0) then
               UL = 1.D0 - mag4/(3.D0*mag2**2)
            else
               UL = 0.D0
            end if
            if (str2 > 0.D0) then
               US = 2.D0 - str4/(str2**2)
            else
               US = 0.D0
            end if
            if (strx2 > 0.D0) then
               USX = 1.D0 - strx4/(3.D0*strx2**2)
            else
               USX = 0.D0
            end if
            if (stry2 > 0.D0) then
               USY = 1.D0 - stry4/(3.D0*stry2**2)
            else
               USY = 0.D0
            end if
            if (nbar2 > 0.D0) then
               UBN = 1.D0 - nbar4/(3.D0*nbar2**2)
            else
               UBN = 0.D0
            end if
            sph = (en2 - en**2)*(L3 - nimp)*beta_ladder(slot)**2

            confmag(slot) = confmag(slot) + absmag
            conf2mag(slot) = conf2mag(slot) + absmag**2
            confmag2(slot) = confmag2(slot) + mag2
            confmag4(slot) = confmag4(slot) + mag4
            if (absmag > 0.D0) then
               conflogmag(slot) = conflogmag(slot) + log(absmag)
            else
               conflogmag(slot) = conflogmag(slot) + log(tiny(1.D0))
            end if
            confsusc(slot) = confsusc(slot) + suscabs
            confUL(slot) = confUL(slot) + UL
            confUS(slot) = confUS(slot) + US
            conf2susc(slot) = conf2susc(slot) + suscabs**2
            conf2UL(slot) = conf2UL(slot) + UL**2
            conf2US(slot) = conf2US(slot) + US**2
            confen(slot) = confen(slot) + en
            confsph(slot) = confsph(slot) + sph
            conf2en(slot) = conf2en(slot) + en**2
            conf2sph(slot) = conf2sph(slot) + sph**2
            confstr2(slot) = confstr2(slot) + str2
            confstr4(slot) = confstr4(slot) + str4
            confabsstr(slot) = confabsstr(slot) + absstr
            confstrx(slot) = confstrx(slot) + totalstrx
            cfg_strx(slot) = totalstrx
            confstrx2(slot) = confstrx2(slot) + strx2
            confstrx4(slot) = confstrx4(slot) + strx4
            confUSX(slot) = confUSX(slot) + USX
            conf2USX(slot) = conf2USX(slot) + USX**2
            confstry(slot) = confstry(slot) + totalstry
            cfg_stry(slot) = totalstry
            cfg_absstr_diag(slot) = absstr
            cfg_stripediff(slot) = stripediff(slot)
            cfg_besten(slot) = best_en(slot)
            cfg_finalen(slot) = final_en(slot)
            confstry2(slot) = confstry2(slot) + stry2
            confstry4(slot) = confstry4(slot) + stry4
            confUSY(slot) = confUSY(slot) + USY
            conf2USY(slot) = conf2USY(slot) + USY**2
            conf2str(slot) = conf2str(slot) + absstr**2
            confsuscstr(slot) = confsuscstr(slot) + susstrip
            conf2suscstr(slot) = conf2suscstr(slot) + susstrip**2
            confnem(slot) = confnem(slot) + nem
            conf2nem(slot) = conf2nem(slot) + nem**2
            confnbar(slot) = confnbar(slot) + nbar
            conf2nbar(slot) = conf2nbar(slot) + nbar**2
            confnbar2(slot) = confnbar2(slot) + nbar2
            confnbar4(slot) = confnbar4(slot) + nbar4
            confUBN(slot) = confUBN(slot) + UBN
            conf2UBN(slot) = conf2UBN(slot) + UBN**2
            cfg_nbar(slot) = nbar
            cfg_nbar2(slot) = nbar2
            cfg_nbar4(slot) = nbar4
            cfg_nbind(slot) = UBN

            if (WRITEHIST .eqv. .true. .and. NMESS > 0) then
               confhist(slot, :) = confhist(slot, :)/NMESS
            end if
            cfg_numflips(slot, :) = confnumflips(slot, :)
         end do
      end subroutine Pt_measurement_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      subroutine Grouped_pt_run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (.not. pt_group_active) return
         call Grouped_initialize_pt_state
         call Grouped_pt_equilibration_loop
         call Grouped_pt_measurement_loop
      end subroutine Grouped_pt_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_initialize_pt_state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: rep
         integer(i4b) :: replica_seed

         cfg_ptswap_attempts(:) = 0
         cfg_ptswap_accepts(:) = 0
         cfg_ptvisitboth = 0
         confhist(:, :) = 0.D0
         hist_en_min_t(:) = 0.D0
         hist_en_max_t(:) = 0.D0
         hist_bin_size_t(:) = 0.D0
         confnumflips(:, :) = 0
         call Mark_vacancy_flipcounts

         if (pt_group_master) then
            do rep = 1, NREPLICA
               rep_at_temp(rep) = rep
               temp_of_rep(rep) = rep
            end do
         end if
         call MPI_BCAST(rep_at_temp, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)
         call MPI_BCAST(temp_of_rep, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)

         rep = pt_group_rank + 1
         replica_seed = IRINIT + (iconf - 1)*NREPLICA + rep
         call kissinit(replica_seed)
         call Initialize_spin_array(s1_rep(0:L3 - 1, rep))
         call kissinit(IRINIT + 100000 + (iconf - 1)*NREPLICA + rep)
      end subroutine Grouped_initialize_pt_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_pt_equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, swap_phase, rep, current_slot, progress_every
         integer(kind=8) :: phase_clock_start
         real(r8b)    :: rep_mag, rep_en, rep_strx, rep_stry

         call Timer_start(phase_clock_start)
         rep = pt_group_rank + 1
         current_slot = temp_of_rep(rep)
         call Initialize_observables(s1_rep(0:L3 - 1, rep), rep_mag, rep_en, rep_strx, rep_stry)
         call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag, rep_en, rep_strx, rep_stry, 'GPTINIT', current_slot, rep)

         swap_phase = 1
         progress_every = max(1, NEQ/10)
         if (pt_group_master) call Write_progress('PT_EQUIL', 0, NEQ)
         do isweep = 1, NEQ
            call Grouped_pt_cycle(swap_phase, rep, current_slot, rep_mag, rep_en, rep_strx, rep_stry)
            if (swap_phase == 1) then
               swap_phase = 2
            else
               swap_phase = 1
            end if
            if (pt_group_master .and. (mod(isweep, progress_every) == 0 .or. isweep == NEQ)) &
     &         call Write_progress('PT_EQUIL', isweep, NEQ)
         end do
         call Write_timing('PT_EQUIL', phase_clock_start)
      end subroutine Grouped_pt_equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_pt_cycle(swap_phase, rep, current_slot, rep_mag, rep_en, rep_strx, rep_stry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in)    :: swap_phase, rep
         integer(i4b), intent(inout) :: current_slot
         real(r8b), intent(inout)    :: rep_mag, rep_en, rep_strx, rep_stry
         integer(i4b)                :: isweep

         do isweep = 1, PTSWAPEVER
            call Metro_sweep(s1_rep(0:L3 - 1, rep), beta_ladder(current_slot), confnumflips(current_slot, :), &
     &         rep_mag, rep_en, rep_strx, rep_stry)
         end do
         call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag, rep_en, rep_strx, rep_stry, 'GPTCYCLE', current_slot, rep)
         call Grouped_pt_attempt_swaps(swap_phase, rep, current_slot, rep_en)
      end subroutine Grouped_pt_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_pt_attempt_swaps(swap_phase, rep, current_slot, rep_en)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in)    :: swap_phase, rep
         integer(i4b), intent(inout) :: current_slot
         real(r8b), intent(in)       :: rep_en
         integer(i4b)                :: slot, slot_start, rep_lo, rep_hi, tmp
         real(r8b)                   :: all_en(NREPLICA), delta, total_en_lo, total_en_hi

         all_en(:) = 0.D0
         call MPI_GATHER(rep_en, 1, MPI_DOUBLE_PRECISION, all_en, 1, MPI_DOUBLE_PRECISION, 0, pt_comm, ierr)

         if (pt_group_master) then
            if (swap_phase == 1) then
               slot_start = 1
            else
               slot_start = 2
            end if
            do slot = slot_start, NREPLICA - 1, 2
               rep_lo = rep_at_temp(slot)
               rep_hi = rep_at_temp(slot + 1)
               cfg_ptswap_attempts(slot) = cfg_ptswap_attempts(slot) + 1

               total_en_lo = all_en(rep_lo)*(L3 - nimp)
               total_en_hi = all_en(rep_hi)*(L3 - nimp)
               delta = (beta_ladder(slot) - beta_ladder(slot + 1))*(total_en_lo - total_en_hi)
               if (delta >= 0.D0 .or. rkiss05() < exp(delta)) then
                  cfg_ptswap_accepts(slot) = cfg_ptswap_accepts(slot) + 1
                  tmp = rep_at_temp(slot)
                  rep_at_temp(slot) = rep_at_temp(slot + 1)
                  rep_at_temp(slot + 1) = tmp
                  temp_of_rep(rep_lo) = slot + 1
                  temp_of_rep(rep_hi) = slot
               end if
            end do
         end if

         call MPI_BCAST(temp_of_rep, NREPLICA, MPI_INTEGER, 0, pt_comm, ierr)
         current_slot = temp_of_rep(rep)
      end subroutine Grouped_pt_attempt_swaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_pt_measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: isweep, slot, rep, current_slot, nmeasure, swap_phase, progress_every
         integer(kind=8) :: phase_clock_start
         real(r8b)    :: rep_mag, rep_en, rep_strx, rep_stry, rep_nbar
         real(r8b)    :: l_mag_sum(NTEMP), l_absmag_sum(NTEMP), l_mag2_sum(NTEMP), l_mag4_sum(NTEMP)
         real(r8b)    :: l_en_sum(NTEMP), l_en2_sum(NTEMP), l_strx_abs_sum(NTEMP), l_stry_abs_sum(NTEMP)
         real(r8b)    :: l_strx2_sum(NTEMP), l_strx4_sum(NTEMP), l_stry2_sum(NTEMP), l_stry4_sum(NTEMP)
         real(r8b)    :: l_str2_sum(NTEMP), l_str4_sum(NTEMP), l_absstr_sum(NTEMP), l_nem_sum(NTEMP)
         real(r8b)    :: l_nbar_sum(NTEMP), l_nbar2_sum(NTEMP), l_nbar4_sum(NTEMP)
         real(r8b)    :: l_best_en(NTEMP), l_final_en(NTEMP)
         real(r8b)    :: mag_sum(NTEMP), absmag_sum(NTEMP), mag2_sum(NTEMP), mag4_sum(NTEMP)
         real(r8b)    :: en_sum(NTEMP), en2_sum(NTEMP), strx_abs_sum(NTEMP), stry_abs_sum(NTEMP)
         real(r8b)    :: strx2_sum(NTEMP), strx4_sum(NTEMP), stry2_sum(NTEMP), stry4_sum(NTEMP)
         real(r8b)    :: str2_sum(NTEMP), str4_sum(NTEMP), absstr_sum(NTEMP), nem_sum(NTEMP)
         real(r8b)    :: nbar_sum(NTEMP), nbar2_sum(NTEMP), nbar4_sum(NTEMP)
         real(r8b)    :: best_en(NTEMP), final_en(NTEMP), stripediff(NTEMP)
         real(r8b)    :: amag, aen, astripx, astripy, anbar, str2, absstr, nem
         real(r8b)    :: mag, absmag, mag2, mag4, en, en2, totalstrx, totalstry
         real(r8b)    :: strx2, strx4, stry2, stry4, str4, suscabs, susstrip, UL, US, USX, USY, sph
         real(r8b)    :: nbar, nbar2, nbar4, UBN
         real(r8b)    :: local_visitboth(1), group_visitboth(1)
         real(r8b)    :: group_locspin_accum(0:L3 - 1), group_locnem_accum(0:L3 - 1)
         real(r8b)    :: local_locnem_samples(1), group_locnem_samples(1)
         logical(ilog) :: visited_cold, visited_hot

         l_mag_sum(:) = 0.D0
         l_absmag_sum(:) = 0.D0
         l_mag2_sum(:) = 0.D0
         l_mag4_sum(:) = 0.D0
         l_en_sum(:) = 0.D0
         l_en2_sum(:) = 0.D0
         l_strx_abs_sum(:) = 0.D0
         l_stry_abs_sum(:) = 0.D0
         l_strx2_sum(:) = 0.D0
         l_strx4_sum(:) = 0.D0
         l_stry2_sum(:) = 0.D0
         l_stry4_sum(:) = 0.D0
         l_str2_sum(:) = 0.D0
         l_str4_sum(:) = 0.D0
         l_absstr_sum(:) = 0.D0
         l_nem_sum(:) = 0.D0
         l_nbar_sum(:) = 0.D0
         l_nbar2_sum(:) = 0.D0
         l_nbar4_sum(:) = 0.D0
         l_best_en(:) = huge(1.D0)
         l_final_en(:) = 0.D0
         cfg_phi(:) = p(:)
         call Reset_local_nematic_average
         call Timer_start(phase_clock_start)

         rep = pt_group_rank + 1
         current_slot = temp_of_rep(rep)
         visited_cold = .false.
         visited_hot = .false.
         call Initialize_observables(s1_rep(0:L3 - 1, rep), rep_mag, rep_en, rep_strx, rep_stry)
         call Spatial_mean_local_nematic(s1_rep(0:L3 - 1, rep), rep_nbar)
         call Validate_observable_state(s1_rep(0:L3 - 1, rep), rep_mag, rep_en, rep_strx, rep_stry, 'GPTMEASINIT', current_slot, rep)

         if (NMESS <= 0) then
            nmeasure = 1
            if (pt_group_master) call Write_progress('PT_MEAS', 0, 0)
            if (WRITELOCNEM .eqv. .true. .and. current_slot == PTMEASINDEX) &
     &         call Accumulate_local_nematic_average(s1_rep(0:L3 - 1, rep))
            call Grouped_accumulate_measurement(rep, current_slot, rep_mag, rep_en, rep_strx, rep_stry, rep_nbar, &
     &         l_mag_sum, l_absmag_sum, l_mag2_sum, l_mag4_sum, l_en_sum, l_en2_sum, l_strx_abs_sum, &
     &         l_stry_abs_sum, l_strx2_sum, l_strx4_sum, l_stry2_sum, l_stry4_sum, l_str2_sum, l_str4_sum, &
     &         l_absstr_sum, l_nem_sum, l_nbar_sum, l_nbar2_sum, l_nbar4_sum, l_best_en)
         else
            nmeasure = NMESS
            swap_phase = 1
            progress_every = max(1, NMESS/10)
            if (pt_group_master) call Write_progress('PT_MEAS', 0, NMESS)
            do isweep = 1, NMESS
               call Grouped_pt_cycle(swap_phase, rep, current_slot, rep_mag, rep_en, rep_strx, rep_stry)
               if (current_slot == 1) visited_cold = .true.
               if (current_slot == NREPLICA) visited_hot = .true.
              call Spatial_mean_local_nematic(s1_rep(0:L3 - 1, rep), rep_nbar)
               if (WRITELOCNEM .eqv. .true. .and. current_slot == PTMEASINDEX) &
     &            call Accumulate_local_nematic_average(s1_rep(0:L3 - 1, rep))
               call Grouped_accumulate_measurement(rep, current_slot, rep_mag, rep_en, rep_strx, rep_stry, rep_nbar, &
     &            l_mag_sum, l_absmag_sum, l_mag2_sum, l_mag4_sum, l_en_sum, l_en2_sum, l_strx_abs_sum, &
     &            l_stry_abs_sum, l_strx2_sum, l_strx4_sum, l_stry2_sum, l_stry4_sum, l_str2_sum, l_str4_sum, &
     &            l_absstr_sum, l_nem_sum, l_nbar_sum, l_nbar2_sum, l_nbar4_sum, l_best_en)

               if (swap_phase == 1) then
                  swap_phase = 2
               else
                  swap_phase = 1
               end if
               if (pt_group_master .and. (mod(isweep, progress_every) == 0 .or. isweep == NMESS)) &
     &            call Write_progress('PT_MEAS', isweep, NMESS)
            end do
         end if
         call Write_timing('PT_MEAS', phase_clock_start)

         l_final_en(:) = 0.D0
         l_final_en(current_slot) = rep_en
         local_visitboth(1) = merge(1.D0, 0.D0, visited_cold .and. visited_hot)
         local_locnem_samples(1) = real(locnem_samples, r8b)

         call MPI_REDUCE(l_mag_sum, mag_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_absmag_sum, absmag_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_mag2_sum, mag2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_mag4_sum, mag4_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_en_sum, en_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_en2_sum, en2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_strx_abs_sum, strx_abs_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_stry_abs_sum, stry_abs_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_strx2_sum, strx2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_strx4_sum, strx4_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_stry2_sum, stry2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_stry4_sum, stry4_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_str2_sum, str2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_str4_sum, str4_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_absstr_sum, absstr_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_nem_sum, nem_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_nbar_sum, nbar_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_nbar2_sum, nbar2_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_nbar4_sum, nbar4_sum, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(l_best_en, best_en, NTEMP, MPI_DOUBLE_PRECISION, MPI_MIN, 0, pt_comm, ierr)
         call MPI_REDUCE(l_final_en, final_en, NTEMP, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         call MPI_REDUCE(local_visitboth, group_visitboth, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         if (WRITELOCNEM .eqv. .true.) then
            call MPI_REDUCE(locspin_accum, group_locspin_accum, L3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
            call MPI_REDUCE(locnem_accum, group_locnem_accum, L3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
            call MPI_REDUCE(local_locnem_samples, group_locnem_samples, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, pt_comm, ierr)
         end if

         if (pt_group_master) then
            if (WRITELOCNEM .eqv. .true.) then
               locspin_accum(:) = group_locspin_accum(:)
               locnem_accum(:) = group_locnem_accum(:)
               locnem_samples = nint(group_locnem_samples(1))
               itemp = PTMEASINDEX
               T = temp_ladder(PTMEASINDEX)
               beta = beta_ladder(PTMEASINDEX)
            end if
            cfg_ptvisitboth = nint(group_visitboth(1))
            do slot = 1, NTEMP
               mag = abs(mag_sum(slot))/nmeasure
               absmag = absmag_sum(slot)/nmeasure
               mag2 = mag2_sum(slot)/nmeasure
               mag4 = mag4_sum(slot)/nmeasure
               en = en_sum(slot)/nmeasure
               en2 = en2_sum(slot)/nmeasure
               totalstrx = strx_abs_sum(slot)/nmeasure
               totalstry = stry_abs_sum(slot)/nmeasure
               strx2 = strx2_sum(slot)/nmeasure
               strx4 = strx4_sum(slot)/nmeasure
               stry2 = stry2_sum(slot)/nmeasure
               stry4 = stry4_sum(slot)/nmeasure
               str2 = str2_sum(slot)/nmeasure
               str4 = str4_sum(slot)/nmeasure
               absstr = absstr_sum(slot)/nmeasure
               nem = nem_sum(slot)/nmeasure
               nbar = nbar_sum(slot)/nmeasure
               nbar2 = nbar2_sum(slot)/nmeasure
               nbar4 = nbar4_sum(slot)/nmeasure
               stripediff(slot) = totalstrx**2 - totalstry**2
               suscabs = (mag2 - absmag**2)*(L3 - nimp)*beta_ladder(slot)
               susstrip = (str2 - absstr**2)*(L3 - nimp)*beta_ladder(slot)
               if (mag2 > 0.D0) then
                  UL = 1.D0 - mag4/(3.D0*mag2**2)
               else
                  UL = 0.D0
               end if
               if (str2 > 0.D0) then
                  US = 2.D0 - str4/(str2**2)
               else
                  US = 0.D0
               end if
               if (strx2 > 0.D0) then
                  USX = 1.D0 - strx4/(3.D0*strx2**2)
               else
                  USX = 0.D0
               end if
               if (stry2 > 0.D0) then
                  USY = 1.D0 - stry4/(3.D0*stry2**2)
               else
                  USY = 0.D0
               end if
               if (nbar2 > 0.D0) then
                  UBN = 1.D0 - nbar4/(3.D0*nbar2**2)
               else
                  UBN = 0.D0
               end if
               sph = (en2 - en**2)*(L3 - nimp)*beta_ladder(slot)**2

               confmag(slot) = confmag(slot) + absmag
               conf2mag(slot) = conf2mag(slot) + absmag**2
               confmag2(slot) = confmag2(slot) + mag2
               confmag4(slot) = confmag4(slot) + mag4
               if (absmag > 0.D0) then
                  conflogmag(slot) = conflogmag(slot) + log(absmag)
               else
                  conflogmag(slot) = conflogmag(slot) + log(tiny(1.D0))
               end if
               confsusc(slot) = confsusc(slot) + suscabs
               confUL(slot) = confUL(slot) + UL
               confUS(slot) = confUS(slot) + US
               conf2susc(slot) = conf2susc(slot) + suscabs**2
               conf2UL(slot) = conf2UL(slot) + UL**2
               conf2US(slot) = conf2US(slot) + US**2
               confen(slot) = confen(slot) + en
               confsph(slot) = confsph(slot) + sph
               conf2en(slot) = conf2en(slot) + en**2
               conf2sph(slot) = conf2sph(slot) + sph**2
               confstr2(slot) = confstr2(slot) + str2
               confstr4(slot) = confstr4(slot) + str4
               confabsstr(slot) = confabsstr(slot) + absstr
               confstrx(slot) = confstrx(slot) + totalstrx
               cfg_strx(slot) = totalstrx
               confstrx2(slot) = confstrx2(slot) + strx2
               confstrx4(slot) = confstrx4(slot) + strx4
               confUSX(slot) = confUSX(slot) + USX
               conf2USX(slot) = conf2USX(slot) + USX**2
               confstry(slot) = confstry(slot) + totalstry
               cfg_stry(slot) = totalstry
               cfg_absstr_diag(slot) = absstr
               cfg_stripediff(slot) = stripediff(slot)
               cfg_besten(slot) = best_en(slot)
               cfg_finalen(slot) = final_en(slot)
               confstry2(slot) = confstry2(slot) + stry2
               confstry4(slot) = confstry4(slot) + stry4
               confUSY(slot) = confUSY(slot) + USY
               conf2USY(slot) = conf2USY(slot) + USY**2
               conf2str(slot) = conf2str(slot) + absstr**2
               confsuscstr(slot) = confsuscstr(slot) + susstrip
               conf2suscstr(slot) = conf2suscstr(slot) + susstrip**2
               confnem(slot) = confnem(slot) + nem
               conf2nem(slot) = conf2nem(slot) + nem**2
               confnbar(slot) = confnbar(slot) + nbar
               conf2nbar(slot) = conf2nbar(slot) + nbar**2
               confnbar2(slot) = confnbar2(slot) + nbar2
               confnbar4(slot) = confnbar4(slot) + nbar4
               confUBN(slot) = confUBN(slot) + UBN
               conf2UBN(slot) = conf2UBN(slot) + UBN**2
               cfg_nbar(slot) = nbar
               cfg_nbar2(slot) = nbar2
               cfg_nbar4(slot) = nbar4
               cfg_nbind(slot) = UBN
            end do
         end if
      end subroutine Grouped_pt_measurement_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Grouped_accumulate_measurement(rep, slot, rep_mag, rep_en, rep_strx, rep_stry, rep_nbar, &
     &   mag_sum, absmag_sum, mag2_sum, mag4_sum, en_sum, en2_sum, strx_abs_sum, stry_abs_sum, &
     &   strx2_sum, strx4_sum, stry2_sum, stry4_sum, str2_sum, str4_sum, absstr_sum, nem_sum, &
     &   nbar_sum, nbar2_sum, nbar4_sum, best_en)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: rep, slot
         real(r8b), intent(in)    :: rep_mag, rep_en, rep_strx, rep_stry, rep_nbar
         real(r8b), intent(inout) :: mag_sum(NTEMP), absmag_sum(NTEMP), mag2_sum(NTEMP), mag4_sum(NTEMP)
         real(r8b), intent(inout) :: en_sum(NTEMP), en2_sum(NTEMP), strx_abs_sum(NTEMP), stry_abs_sum(NTEMP)
         real(r8b), intent(inout) :: strx2_sum(NTEMP), strx4_sum(NTEMP), stry2_sum(NTEMP), stry4_sum(NTEMP)
         real(r8b), intent(inout) :: str2_sum(NTEMP), str4_sum(NTEMP), absstr_sum(NTEMP), nem_sum(NTEMP)
         real(r8b), intent(inout) :: nbar_sum(NTEMP), nbar2_sum(NTEMP), nbar4_sum(NTEMP), best_en(NTEMP)
         real(r8b)                :: str2, absstr, nem

         str2 = rep_strx**2 + rep_stry**2
         absstr = sqrt(str2)
         nem = abs(rep_stry**2 - rep_strx**2)
         mag_sum(slot) = mag_sum(slot) + rep_mag
         absmag_sum(slot) = absmag_sum(slot) + abs(rep_mag)
         mag2_sum(slot) = mag2_sum(slot) + rep_mag**2
         mag4_sum(slot) = mag4_sum(slot) + rep_mag**4
         en_sum(slot) = en_sum(slot) + rep_en
         en2_sum(slot) = en2_sum(slot) + rep_en**2
         strx_abs_sum(slot) = strx_abs_sum(slot) + abs(rep_strx)
         stry_abs_sum(slot) = stry_abs_sum(slot) + abs(rep_stry)
         strx2_sum(slot) = strx2_sum(slot) + rep_strx**2
         strx4_sum(slot) = strx4_sum(slot) + rep_strx**4
         stry2_sum(slot) = stry2_sum(slot) + rep_stry**2
         stry4_sum(slot) = stry4_sum(slot) + rep_stry**4
         str2_sum(slot) = str2_sum(slot) + str2
         str4_sum(slot) = str4_sum(slot) + str2**2
         absstr_sum(slot) = absstr_sum(slot) + absstr
         nem_sum(slot) = nem_sum(slot) + nem
         nbar_sum(slot) = nbar_sum(slot) + rep_nbar
         nbar2_sum(slot) = nbar2_sum(slot) + rep_nbar**2
         nbar4_sum(slot) = nbar4_sum(slot) + rep_nbar**4
         best_en(slot) = min(best_en(slot), rep_en)
      end subroutine Grouped_accumulate_measurement
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Snapshot_pt_target_replica
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: rep

         rep = rep_at_temp(PTMEASINDEX)
         s1(:) = s1_rep(:, rep)
         itemp = PTMEASINDEX
         T = temp_ladder(PTMEASINDEX)
         beta = beta_ladder(PTMEASINDEX)
      end subroutine Snapshot_pt_target_replica

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Equilibration_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         real(r8b)       :: curmag, curen, curstrx, curstry

         call Initialize_observables(s1, curmag, curen, curstrx, curstry)
         call Validate_observable_state(s1, curmag, curen, curstrx, curstry, 'LOCALEQINIT', 0, itemp)

         do isweep = 1, NEQ
            call Metro_sweep(s1, beta, confnumflips(itemp, :), curmag, curen, curstrx, curstry)
            call Validate_observable_state(s1, curmag, curen, curstrx, curstry, 'LOCALEQ', isweep, itemp)
            ! corner_sweep is intentionally left disabled here: as written below it
            ! does not include the active phi coupling and is therefore not used as
            ! a baseline-equilibration move under the current Hamiltonian.
         end do

      end subroutine Equilibration_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Measurement_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         real(r8b)       :: curmag, curen         ! magnetization and energy for a single sweep
         real(r8b)       :: curstrx, curstry         ! stripex and stripey for a single sweep
         real(r8b)       :: curnbar
         real(r8b)       :: sumstrx2, sumstrx4, sumstry2, sumstry4
         real(r8b)       :: avgstrx2, avgstrx4, avgstry2, avgstry4
         real(r8b)       :: bindx, bindy
         logical(ilog)   :: do_bind_time
         character(len=21) :: bt_trace_file

         call Reset_local_nematic_average
         call Initialize_observables(s1, curmag, curen, curstrx, curstry)
         call Validate_observable_state(s1, curmag, curen, curstrx, curstry, 'LOCALINIT', 0, itemp)

         if (NMESS <= 0) then
            call Accumulate_local_nematic_average(s1)
            return
         end if

         do_bind_time = WRITEBINDTIME .and. (iconf <= NTRACECONF) .and. (abs(T - TMIN) < 1.0D-8)
         if (do_bind_time) then
            write (bt_trace_file, '("cfbt",I4.4,I4.4,"_c",I3.3,".dat")') L, H, iconf
            open (77, file=bt_trace_file, status='replace')
            write (77, *) 'iconf T'
            write (77, '(I4,1x,E12.6)') iconf, T
            write (77, *) 'isweep, binderx, bindery'
            sumstrx2 = 0.D0
            sumstrx4 = 0.D0
            sumstry2 = 0.D0
            sumstry4 = 0.D0
         end if

         do isweep = 1, NMESS
            call Metro_sweep(s1, beta, confnumflips(itemp, :), curmag, curen, curstrx, curstry)
            call Validate_observable_state(s1, curmag, curen, curstrx, curstry, 'LOCALMEAS', isweep, itemp)
            call Accumulate_local_nematic_average(s1)
            call Spatial_mean_local_nematic(s1, curnbar)
            maglist(isweep) = curmag
            enlist(isweep) = curen
            strxlist(isweep) = curstrx
            strylist(isweep) = curstry
            nbarlist(isweep) = curnbar

            if (do_bind_time) then
               sumstrx2 = sumstrx2 + curstrx**2
               sumstrx4 = sumstrx4 + curstrx**4
               sumstry2 = sumstry2 + curstry**2
               sumstry4 = sumstry4 + curstry**4

               avgstrx2 = sumstrx2/isweep
               avgstrx4 = sumstrx4/isweep
               avgstry2 = sumstry2/isweep
               avgstry4 = sumstry4/isweep

               if (avgstrx2 > 0.D0) then
                  bindx = 1.D0 - avgstrx4/(3.0D0*avgstrx2**2)
               else
                  bindx = 0.D0
               end if
               if (avgstry2 > 0.D0) then
                  bindy = 1.D0 - avgstry4/(3.0D0*avgstry2**2)
               else
                  bindy = 0.D0
               end if

               write (77, '(I8,1x,E12.6,1x,E12.6)') isweep, bindx, bindy
            end if

         end do
         if (do_bind_time) then
            close (77)
         end if

      end subroutine Measurement_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Process_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)    :: isweep                   ! sweep counter
         real(r8b)       :: amag, aen                ! current magnetazation, energy
         real(r8b)       :: mag, absmag, mag2, mag4     ! average magnetization, abs. value, square , 4th power
         real(r8b)       :: astripx, astripy, totalstrx, totalstry     ! strip x, strip y
         real(r8b)       :: strx2, stry2, strx4, stry4  ! strip x^2, strip y^2, strip x^4, strip y^4
         real(r8b)       :: str2, str4, absstr     ! average strip square, 4th power        ,abs.value
         real(r8b)       :: susccl, suscabs, susstrip          ! susceptibilities calculated from clustersize and mag2 and strip
         real(r8b)       :: UL, US                       ! Binder parameter
         real(r8b)       :: USX, USY                     ! Binder parameter for strip x and y
         real(r8b)       :: en, en2, sph             ! average energy, square, specific heat
         real(r8b)       :: nem                      ! nematic OP
         real(r8b)       :: nbar, nbar2, nbar4, UBN  ! spatial mean local nematic and Binder
         real(r8b)       :: hrange                   ! histogram range padding
         real(r8b)       :: measured_mag, measured_en, measured_strx, measured_stry, measured_nbar
         integer(i4b)    :: histindex                ! index for histogram
         real(r8b)       :: besten, finalen, stripediff, anbar

         mag = 0.D0
         absmag = 0.D0
         mag2 = 0.D0
         mag4 = 0.D0
         en = 0.D0
         en2 = 0.D0
         astripx = 0.D0
         totalstrx = 0.D0
         strx2 = 0.D0
         strx4 = 0.D0
         astripy = 0.D0
         totalstry = 0.D0
         stry2 = 0.D0
         stry4 = 0.D0
         str2 = 0.D0
         str4 = 0.D0
         absstr = 0.D0
         nem = 0.D0
         nbar = 0.D0
         nbar2 = 0.D0
         nbar4 = 0.D0
         UBN = 0.D0
         besten = 0.D0
         finalen = 0.D0
         stripediff = 0.D0

         if (NMESS <= 0) then
            call Measurement(s1, measured_mag, measured_en, measured_strx, measured_stry)
            call Spatial_mean_local_nematic(s1, measured_nbar)
            mag = abs(measured_mag)
            absmag = abs(measured_mag)
            mag2 = measured_mag**2
            mag4 = measured_mag**4
            en = measured_en
            en2 = measured_en**2
            totalstrx = abs(measured_strx)
            totalstry = abs(measured_stry)
            strx2 = measured_strx**2
            strx4 = measured_strx**4
            stry2 = measured_stry**2
            stry4 = measured_stry**4
            str2 = measured_strx**2 + measured_stry**2
            str4 = (measured_strx**2 + measured_stry**2)**2
            absstr = sqrt(measured_strx**2 + measured_stry**2)
            nem = abs(measured_stry**2 - measured_strx**2)
            nbar = measured_nbar
            nbar2 = measured_nbar**2
            nbar4 = measured_nbar**4
            besten = measured_en
            finalen = measured_en
            stripediff = totalstrx**2 - totalstry**2
            if (WRITEHIST .eqv. .true.) then
               confhist(itemp, :) = 0.D0
            end if
         else
            ! initialize histogram for this temperature and disorder config
            ! set bin ranges once per temperature (first time we see this itemp)
            if (WRITEHIST .eqv. .true.) then
               confhist(itemp, :) = 0.D0
               if (hist_bin_size_t(itemp) <= 0.D0) then
                  hist_en_min_t(itemp) = minval(enlist)
                  hist_en_max_t(itemp) = maxval(enlist)
                  if (hist_en_max_t(itemp) > hist_en_min_t(itemp)) then
                     hrange = hist_en_max_t(itemp) - hist_en_min_t(itemp)
                     hist_en_min_t(itemp) = hist_en_min_t(itemp) - 0.50D0*hrange
                     hist_en_max_t(itemp) = hist_en_max_t(itemp) + 0.50D0*hrange
                  end if
                  if (hist_en_max_t(itemp) <= hist_en_min_t(itemp)) then
                     hist_en_min_t(itemp) = hist_en_min_t(itemp) - 1.0D-6
                     hist_en_max_t(itemp) = hist_en_max_t(itemp) + 1.0D-6
                  end if
                  hist_bin_size_t(itemp) = (hist_en_max_t(itemp) - hist_en_min_t(itemp))/BINS
               end if
            end if

            besten = minval(enlist)
            finalen = enlist(NMESS)

            do isweep = 1, NMESS
               amag = maglist(isweep)
               mag = mag + amag
               astripx = strxlist(isweep)
               totalstrx = totalstrx + abs(astripx)
               strx2 = strx2 + astripx**2
               strx4 = strx4 + astripx**4
               astripy = strylist(isweep)
               totalstry = totalstry + abs(astripy)
               stry2 = stry2 + astripy**2
               stry4 = stry4 + astripy**4
               str2 = str2 + astripx**2 + astripy**2
               str4 = str4 + (astripx**2 + astripy**2)**2
               absstr = absstr + sqrt(astripx**2 + astripy**2)
               absmag = absmag + abs(amag)
               mag2 = mag2 + amag**2
               mag4 = mag4 + amag**4
               nem = nem + abs(astripy**2 - astripx**2)
               anbar = nbarlist(isweep)
               nbar = nbar + anbar
               nbar2 = nbar2 + anbar**2
               nbar4 = nbar4 + anbar**4

               aen = enlist(isweep)
               en = en + aen
               en2 = en2 + aen**2

               if (WRITEHIST .eqv. .true.) then
                  histindex = NINT((aen - hist_en_min_t(itemp))/hist_bin_size_t(itemp))
                  if (histindex .ge. 1 .and. histindex .le. BINS) then
                     confhist(itemp, histindex) = confhist(itemp, histindex) + 1
                  end if
               end if

            end do

            mag = abs(mag)/NMESS
            absmag = absmag/NMESS
            mag2 = mag2/NMESS
            mag4 = mag4/NMESS
            en = en/NMESS
            en2 = en2/NMESS
            str2 = str2/NMESS
            str4 = str4/NMESS
            absstr = absstr/NMESS
            nem = nem/NMESS
            nbar = nbar/NMESS
            nbar2 = nbar2/NMESS
            nbar4 = nbar4/NMESS
            totalstrx = totalstrx/NMESS
            totalstry = totalstry/NMESS
            strx2 = strx2/NMESS
            strx4 = strx4/NMESS
            stry2 = stry2/NMESS
            stry4 = stry4/NMESS
            stripediff = totalstrx**2 - totalstry**2
         end if

!          print *, 'nem=',nem

         suscabs = (mag2 - absmag**2)*(L3 - nimp)*beta
         susccl = (avclsize/(L3) - absmag**2)*(L3 - nimp)*beta
         susstrip = (str2 - absstr**2)*(L3 - nimp)*beta
         
         ! zero gaurds
         if (mag2 > 0.D0) then
            UL = 1.D0 - mag4/(3.D0*mag2**2)
         else
            UL = 0.D0
         end if
         if (str2 > 0.D0) then
            US = 2.D0 - str4/(str2**2)
         else
            US = 0.D0
         end if
         if (strx2 > 0.D0) then
            USX = 1.D0 - strx4/(3.0D0*strx2**2)
         else
            USX = 0.D0
         end if
         if (stry2 > 0.D0) then
            USY = 1.D0 - stry4/(3.0D0*stry2**2)
         else
            USY = 0.D0
         end if
         if (nbar2 > 0.D0) then
            UBN = 1.D0 - nbar4/(3.D0*nbar2**2)
         else
            UBN = 0.D0
         end if
         sph = (en2 - en**2)*(L3 - nimp)*beta**2

         confmag(itemp) = confmag(itemp) + absmag
         conf2mag(itemp) = conf2mag(itemp) + absmag**2
         confmag2(itemp) = confmag2(itemp) + mag2
         confmag4(itemp) = confmag4(itemp) + mag4
         if (absmag > 0.D0) then
            conflogmag(itemp) = conflogmag(itemp) + log(absmag)
         else
            conflogmag(itemp) = conflogmag(itemp) + log(tiny(1.D0))
         end if
         confsusc(itemp) = confsusc(itemp) + suscabs
         confUL(itemp) = confUL(itemp) + UL
         confUS(itemp) = confUS(itemp) + US
         conf2susc(itemp) = conf2susc(itemp) + suscabs**2
         conf2UL(itemp) = conf2UL(itemp) + UL**2
         conf2US(itemp) = conf2US(itemp) + US**2
         confen(itemp) = confen(itemp) + en
         confsph(itemp) = confsph(itemp) + sph
         conf2en(itemp) = conf2en(itemp) + en**2
         conf2sph(itemp) = conf2sph(itemp) + sph**2
         confstr2(itemp) = confstr2(itemp) + str2
         confstr4(itemp) = confstr4(itemp) + str4
         confabsstr(itemp) = confabsstr(itemp) + absstr
         confstrx(itemp) = confstrx(itemp) + totalstrx
         cfg_strx(itemp) = totalstrx
         confstrx2(itemp) = confstrx2(itemp) + strx2
         confstrx4(itemp) = confstrx4(itemp) + strx4
         confUSX(itemp) = confUSX(itemp) + USX
         conf2USX(itemp) = conf2USX(itemp) + USX**2
         confstry(itemp) = confstry(itemp) + totalstry
         cfg_stry(itemp) = totalstry
         cfg_absstr_diag(itemp) = absstr
         cfg_stripediff(itemp) = stripediff
         cfg_besten(itemp) = besten
         cfg_finalen(itemp) = finalen
         confstry2(itemp) = confstry2(itemp) + stry2
         confstry4(itemp) = confstry4(itemp) + stry4
         confUSY(itemp) = confUSY(itemp) + USY
         conf2USY(itemp) = conf2USY(itemp) + USY**2
         conf2str(itemp) = conf2str(itemp) + absstr**2
         confsuscstr(itemp) = confsuscstr(itemp) + susstrip
         conf2suscstr(itemp) = conf2suscstr(itemp) + susstrip**2
         confnem(itemp) = confnem(itemp) + nem
         conf2nem(itemp) = conf2nem(itemp) + nem**2
         confnbar(itemp) = confnbar(itemp) + nbar
         conf2nbar(itemp) = conf2nbar(itemp) + nbar**2
         confnbar2(itemp) = confnbar2(itemp) + nbar2
         confnbar4(itemp) = confnbar4(itemp) + nbar4
         confUBN(itemp) = confUBN(itemp) + UBN
         conf2UBN(itemp) = conf2UBN(itemp) + UBN**2
         cfg_nbar(itemp) = nbar
         cfg_nbar2(itemp) = nbar2
         cfg_nbar4(itemp) = nbar4
         cfg_nbind(itemp) = UBN

         if (WRITEHIST .eqv. .true. .and. NMESS > 0) then
            confhist(itemp, :) = confhist(itemp, :)/NMESS
         end if
         cfg_numflips(itemp, :) = confnumflips(itemp, :)
         cfg_phi(:) = p(:)

      end subroutine Process_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Build_sweep_site_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b) :: idx, jsite, tmp

         do idx = 0, L3 - 1
            site_order(idx) = idx
         end do

         if (UPDATEORDER == 'RND') then
            ! Generate one fresh Fisher-Yates permutation per sweep from the
            ! existing simulation RNG so randomized sweeps stay seed-reproducible.
            do idx = L3 - 1, 1, -1
               jsite = int((idx + 1)*rkiss05())
               tmp = site_order(idx)
               site_order(idx) = site_order(jsite)
               site_order(jsite) = tmp
            end do
         end if
      end subroutine Build_sweep_site_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Compute_local_field(spins, site, field, n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Local field for the active Hamiltonian, including the phi coupling.
!!    Metro_sweep uses this helper so the sequential and randomized sweep
!!    variants share identical acceptance physics.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in)  :: spins(0:L3 - 1)
         integer(i4b), intent(in) :: site
         real(r8b), intent(out)   :: field
         real(r8b), intent(out)   :: n

         field = 0.0D0
         n = 0.0D0
         if (occu(m1(site)) == 1) then
            field = field + spins(m1(site))*J1v*jv(site)
            n = n - spins(m1(site))
            field = field - p(m1(site))*spins(m1(site))
         end if
         if (occu(m2(site)) == 1) then
            field = field + spins(m2(site))*J1v*jv(m2(site))
            n = n - spins(m2(site))
            field = field - p(m2(site))*spins(m2(site))
         end if
         if (occu(m3(site)) == 1) then
            field = field + spins(m3(site))*J1h*jh(site)
            n = n + spins(m3(site))
            field = field + p(m3(site))*spins(m3(site))
         end if
         if (occu(m4(site)) == 1) then
            field = field + spins(m4(site))*J1h*jh(m4(site))
            n = n + spins(m4(site))
            field = field + p(m4(site))*spins(m4(site))
         end if

         if (occu(m5(site)) == 1) field = field + spins(m5(site))*J2
         if (occu(m6(site)) == 1) field = field + spins(m6(site))*J2
         if (occu(m7(site)) == 1) field = field + spins(m7(site))*J2
         if (occu(m8(site)) == 1) field = field + spins(m8(site))*J2

         if (occu(m9(site)) == 1) field = field + spins(m9(site))*J3
         if (occu(m10(site)) == 1) field = field + spins(m10(site))*J3

         field = field + (p(site)*n)
      end subroutine Compute_local_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Metro_sweep(spins, beta_local, flipcounts, curmag, curen, curstrx, curstry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Performs a Metropolis sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(inout) :: spins(0:L3 - 1)
         real(r8b), intent(in)       :: beta_local
         integer(i4b), intent(inout) :: flipcounts(0:L3 - 1)
         real(r8b), intent(inout)    :: curmag, curen, curstrx, curstry
         real(r8b)   :: field
         real(r8b)   :: n     ! nematic
         real(r8b)   :: norm, delta_mag, delta_en, delta_strx, delta_stry
         integer(i4b) :: iorder
         integer(i1b)  :: sold

         call Build_sweep_site_order
         norm = real(L3 - nimp, r8b)

         do iorder = 0, L3 - 1
         isite = site_order(iorder)
         if (occu(isite) == 1) then
            ! Evaluate the current site's local energy change.
            call Compute_local_field(spins, isite, field, n)
            sold = spins(isite)
            delta_mag = -2.D0*real(sold, r8b)/norm
            delta_en = 2.D0*real(sold, r8b)*field/norm
            delta_strx = delta_mag*stripe_parity_x(isite)
            delta_stry = delta_mag*stripe_parity_y(isite)

            if (spins(isite) == 1) then
               if (field .le. 0) then
                  spins(isite) = -spins(isite)
                  flipcounts(isite) = flipcounts(isite) + 1
                  curmag = curmag + delta_mag
                  curen = curen + delta_en
                  curstrx = curstrx + delta_strx
                  curstry = curstry + delta_stry
               elseif (rkiss05() .lt. exp(-2.D0*field*beta_local)) then
                  spins(isite) = -spins(isite)
                  flipcounts(isite) = flipcounts(isite) + 1
                  curmag = curmag + delta_mag
                  curen = curen + delta_en
                  curstrx = curstrx + delta_strx
                  curstry = curstry + delta_stry
               end if
            else
               if (field .ge. 0) then
                  spins(isite) = -spins(isite)
                  flipcounts(isite) = flipcounts(isite) + 1
                  curmag = curmag + delta_mag
                  curen = curen + delta_en
                  curstrx = curstrx + delta_strx
                  curstry = curstry + delta_stry
               elseif (rkiss05() .lt. exp(2.D0*field*beta_local)) then
                  spins(isite) = -spins(isite)
                  flipcounts(isite) = flipcounts(isite) + 1
                  curmag = curmag + delta_mag
                  curen = curen + delta_en
                  curstrx = curstrx + delta_strx
                  curstry = curstry + delta_stry
               end if
            end if
         end if
         end do

      end subroutine Metro_sweep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine corner_sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Performs a corner Metropolis sweep
!!    Left inactive: this move does not currently include the phi term used
!!    by Metro_sweep, so enabling it would sample a different Hamiltonian.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b)   :: field

         do isite = 0, L3 - 1                                       ! diagonal exchange top right - bottom left
         if (occu(isite) == 1) then
            if (occu(m5(isite)) == 1) then
               if (s1(isite)*s1(m5(isite)) == -1) then
                  field = 0.0
                  if (occu(m1(isite)) == 1) then
                     field = field + s1(m1(isite))*J1v*jv(isite)
                  end if
                  if (occu(m2(isite)) == 1) then
                     field = field + s1(m2(isite))*J1v*jv(m2(isite))
                  end if
                  if (occu(m3(isite)) == 1) then
                     field = field + s1(m3(isite))*J1h*jh(isite)
                  end if
                  if (occu(m4(isite)) == 1) then
                     field = field + s1(m4(isite))*J1h*jh(m4(isite))
                  end if

                  !         if (occu(m5(isite))==1) then
                  !            field=field+s1(m5(isite))*J2
                  !         endif
                  if (occu(m6(isite)) == 1) then
                     field = field + s1(m6(isite))*J2
                  end if
                  if (occu(m7(isite)) == 1) then
                     field = field + s1(m7(isite))*J2
                  end if
                  if (occu(m8(isite)) == 1) then
                     field = field + s1(m8(isite))*J2
                  end if

                  if (occu(m9(isite)) == 1) then
                     field = field + s1(m9(isite))*J3
                  end if
                  if (occu(m10(isite)) == 1) then
                     field = field + s1(m10(isite))*J3
                  end if

                  if (occu(m1(m5(isite))) == 1) then
                     field = field - s1(m1(m5(isite)))*J1v*jv(m5(isite))
                  end if
                  if (occu(m2(m5(isite))) == 1) then
                     field = field - s1(m2(m5(isite)))*J1v*jv(m2(m5(isite)))
                  end if
                  if (occu(m3(m5(isite))) == 1) then
                     field = field - s1(m3(m5(isite)))*J1h*jh(m5(isite))
                  end if
                  if (occu(m4(m5(isite))) == 1) then
                     field = field - s1(m4(m5(isite)))*J1h*jh(m4(m5(isite)))
                  end if

                  if (occu(m5(m5(isite))) == 1) then
                     field = field - s1(m5(m5(isite)))*J2
                  end if
                  if (occu(m6(m5(isite))) == 1) then
                     field = field - s1(m6(m5(isite)))*J2
                  end if
                  if (occu(m7(m5(isite))) == 1) then
                     field = field - s1(m7(m5(isite)))*J2
                  end if
                  !         if (occu(m8(m5(isite)))==1) then
                  !            field=field-s1(m8(

                  if (occu(m9(m5(isite))) == 1) then
                     field = field - s1(m9(m5(isite)))*J3
                  end if
                  if (occu(m10(m5(isite))) == 1) then
                     field = field - s1(m10(m5(isite)))*J3
                  end if

                  if (s1(isite) == 1) then
                     if (field .le. 0) then
                        s1(isite) = -s1(isite)
                        s1(m5(isite)) = -s1(m5(isite))
                     elseif (rkiss05() .lt. exp(-2.D0*field*beta)) then
                        s1(isite) = -s1(isite)
                        s1(m5(isite)) = -s1(m5(isite))
                     end if
                  else
                     if (field .ge. 0) then
                        s1(isite) = -s1(isite)
                        s1(m5(isite)) = -s1(m5(isite))
                     elseif (rkiss05() .lt. exp(2.D0*field*beta)) then
                        s1(isite) = -s1(isite)
                        s1(m5(isite)) = -s1(m5(isite))
                     end if
                  end if
               end if
            end if
         end if
         end do

         do isite = 0, L3 - 1                                       ! diagonal exchange top left - bottom right
         if (occu(isite) == 1) then
            if (occu(m7(isite)) == 1) then
            if (s1(isite)*s1(m7(isite)) == -1) then
               field = 0.0
               if (occu(m1(isite)) == 1) then
                  field = field + s1(m1(isite))*J1v*jv(isite)
               end if
               if (occu(m2(isite)) == 1) then
                  field = field + s1(m2(isite))*J1v*jv(m2(isite))
               end if
               if (occu(m3(isite)) == 1) then
                  field = field + s1(m3(isite))*J1h*jh(isite)
               end if
               if (occu(m4(isite)) == 1) then
                  field = field + s1(m4(isite))*J1h*jh(m4(isite))
               end if

               if (occu(m5(isite)) == 1) then
                  field = field + s1(m5(isite))*J2
               end if
               if (occu(m6(isite)) == 1) then
                  field = field + s1(m6(isite))*J2
               end if
!         if (occu(m7(isite))==1) then
!            field=field+s1(m7(isite))*J2
!         endif
               if (occu(m8(isite)) == 1) then
                  field = field + s1(m8(isite))*J2
               end if

               if (occu(m9(isite)) == 1) then
                  field = field + s1(m9(isite))*J3
               end if
               if (occu(m10(isite)) == 1) then
                  field = field + s1(m10(isite))*J3
               end if

               if (occu(m1(m7(isite))) == 1) then
                  field = field - s1(m1(m7(isite)))*J1v*jv(m7(isite))
               end if
               if (occu(m2(m7(isite))) == 1) then
                  field = field - s1(m2(m7(isite)))*J1v*jv(m2(m7(isite)))
               end if
               if (occu(m3(m7(isite))) == 1) then
                  field = field - s1(m3(m7(isite)))*J1h*jh(m7(isite))
               end if
               if (occu(m4(m7(isite))) == 1) then
                  field = field - s1(m4(m7(isite)))*J1h*jh(m4(m7(isite)))
               end if

               if (occu(m5(m7(isite))) == 1) then
                  field = field - s1(m5(m7(isite)))*J2
               end if
!         if (occu(m6(m7(isite)))==1) then
!            field=field-s1(m6(m7(isite)))*J2
!         endif
               if (occu(m7(m7(isite))) == 1) then
                  field = field - s1(m7(m7(isite)))*J2
               end if
               if (occu(m8(m7(isite))) == 1) then
                  field = field - s1(m8(m7(isite)))*J2
               end if

               if (occu(m9(m7(isite))) == 1) then
                  field = field - s1(m9(m7(isite)))*J3
               end if
               if (occu(m10(m7(isite))) == 1) then
                  field = field - s1(m10(m7(isite)))*J3
               end if

               if (s1(isite) == 1) then
                  if (field .le. 0) then
                     s1(isite) = -s1(isite)
                     s1(m7(isite)) = -s1(m7(isite))
                  elseif (rkiss05() .lt. exp(-2.D0*field*beta)) then
                     s1(isite) = -s1(isite)
                     s1(m7(isite)) = -s1(m7(isite))
                  end if
               else
                  if (field .ge. 0) then
                     s1(isite) = -s1(isite)
                     s1(m7(isite)) = -s1(m7(isite))
                  elseif (rkiss05() .lt. exp(2.D0*field*beta)) then
                     s1(isite) = -s1(isite)
                     s1(m7(isite)) = -s1(m7(isite))
                  end if
               end if
            end if
            end if
         end if
         end do

      end subroutine corner_sweep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Measurement(spins, curmag, curen, curstrx, curstry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Calculates energy and magnetization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i1b), intent(in)  :: spins(0:L3 - 1)
         real(r8b)       :: imag, ien               ! magnetization and full measured energy numerator
         real(r8b)       :: curmag, curen, curstrx, curstry
         real(r8b)       :: stripx, stripy     ! strip x, strip y
         real(r8b)       :: field, n
         integer(i4b)    :: iz                  ! site index

         ien = 0
         imag = 0
         stripx = 0
         stripy = 0

         do iz = 0, H - 1
         do isite = 0, L2 - 1
            if (occu(isite + L2*iz) == 1) then
               imag = imag + spins(isite + L2*iz)
               call Compute_local_field(spins, isite + L2*iz, field, n)
               ien = ien - spins(isite + L2*iz)*field
               stripx = stripx + spins(isite + L2*iz)*(-1)**mod(isite, L)
               stripy = stripy + spins(isite + L2*iz)*(-1)**(isite/L)
            end if
         end do
         end do

         curmag = (1.D0*imag)/(L3 - nimp)
         ! The local field is the same one used in Metropolis acceptance, so
         ! the 1/2 prefactor here reports the full simulated Hamiltonian.
         curen = (0.5D0*ien)/(L3 - nimp)
         curstrx = (1.D0*stripx)/(L3 - nimp)
         curstry = (1.D0*stripy)/(L3 - nimp)

      end subroutine Measurement

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Local_nematic_OP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Calculates and outputs local nematic order parameters on a cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b)     :: ix, iy, iz                   ! x and y coordinates
         character(len=20) :: filename
         real(r8b)         :: sample_norm, mean_locnem, expected_mean_locnem
         integer(i4b)      :: occ_count

         filename = 'locnem0000.dat'

         write (filename(7:10), '(I4.4)') itemp




         open (8, file=filename, status='replace')
         rewind (8)
         sample_norm = real(max(1, locnem_samples), r8b)
         occ_count = max(1, L3 - nimp)
         mean_locnem = sum(locnem_accum(:))/(sample_norm*real(occ_count, r8b))
         expected_mean_locnem = 4.D0*(cfg_stry(itemp)**2 - cfg_strx(itemp)**2)
         write (8, '(A,E24.16)') '# Temperature = ', T
         write (8, '(A,I0)') '# disorder configuration = ', iconf
         write (8, '(A,2(E24.16,1x))') '# stripe magnitudes = ', cfg_strx(itemp), cfg_stry(itemp)
         write (8, '(A,E24.16)') '# expected mean locnem = ', expected_mean_locnem
         write (8, '(A,E24.16)') '# measured mean locnem = ', mean_locnem
         write (8, '(A)') '# ix iy iz isite occu <s(isite)> <locnem(isite)>'
         do iz = 0, H - 1
         do iy = 0, L - 1
         do ix = 0, L - 1
            isite = L*iy + ix + L2*iz
            write (8, '(4(I6,1x),I3,1x,2(E18.10,1x))') ix, iy, iz, isite, occu(isite), &
               locspin_accum(isite)/sample_norm, locnem_accum(isite)/sample_norm
         end do
         end do
         end do
         close (8)

      end subroutine Local_nematic_OP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Domain_energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Calculates energy gain due to domain formation by flipping every second spin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(r8b)  :: mag0, en0, strx0, stry0
         real(r8b)  :: mag1, en1, strx1, stry1
         integer(i4b)     :: ix, iy                   ! x and y coordinates
         integer(i4b)     :: locnem(0:L2 - 1)

         call measurement(s1, mag0, en0, strx0, stry0)
         do iy = 0, L - 1
         do ix = 0, L - 1
            isite = L*iy + ix
            s1(isite) = s1(isite)*(-1)**ix*(-1)**iy
         end do
         end do
         call measurement(s1, mag1, en1, strx1, stry1)

         print *, 'strx0,stry0,strx1,stry1', strx0, stry0, strx1, stry1
         print *, 'domain energy per site', 0.5D0*(en1 - en0)

         do iy = 0, L - 1                                         ! flip back
            do ix = 0, L - 1
               isite = L*iy + ix
               s1(isite) = s1(isite)*(-1)**ix*(-1)**iy
            end do
         end do

      end subroutine Domain_energy

      subroutine Collect_data(active_config)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         logical(ilog), intent(in) :: active_config
         integer(i4b), parameter :: BASEFIELDS = 47
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
               transdata(1:NTEMP) = confmag(:)
               transdata(NTEMP + 1:2*NTEMP) = confmag2(:)
               transdata(2*NTEMP + 1:3*NTEMP) = confmag4(:)
               transdata(3*NTEMP + 1:4*NTEMP) = confsusc(:)
               transdata(4*NTEMP + 1:5*NTEMP) = conf2susc(:)
               transdata(5*NTEMP + 1:6*NTEMP) = confUL(:)
               transdata(6*NTEMP + 1:7*NTEMP) = conf2UL(:)
               transdata(7*NTEMP + 1:8*NTEMP) = conf2mag(:)
               transdata(8*NTEMP + 1:9*NTEMP) = conflogmag(:)
               transdata(9*NTEMP + 1:10*NTEMP) = confen(:)
               transdata(10*NTEMP + 1:11*NTEMP) = conf2en(:)
               transdata(11*NTEMP + 1:12*NTEMP) = confsph(:)
               transdata(12*NTEMP + 1:13*NTEMP) = conf2sph(:)
               transdata(13*NTEMP + 1:14*NTEMP) = confabsstr(:)
               transdata(14*NTEMP + 1:15*NTEMP) = confstrx(:)
               transdata(15*NTEMP + 1:16*NTEMP) = confstry(:)
               transdata(16*NTEMP + 1:17*NTEMP) = confstr2(:)
               transdata(17*NTEMP + 1:18*NTEMP) = confstr4(:)
               transdata(18*NTEMP + 1:19*NTEMP) = confsuscstr(:)
               transdata(19*NTEMP + 1:20*NTEMP) = conf2suscstr(:)
               transdata(20*NTEMP + 1:21*NTEMP) = confUS(:)
               transdata(21*NTEMP + 1:22*NTEMP) = conf2US(:)
               transdata(22*NTEMP + 1:23*NTEMP) = confnem(:)
               transdata(23*NTEMP + 1:24*NTEMP) = conf2nem(:)
               transdata(24*NTEMP + 1:25*NTEMP) = conf2str(:)
               transdata(25*NTEMP + 1:26*NTEMP) = confstrx2(:)
               transdata(26*NTEMP + 1:27*NTEMP) = confstrx4(:)
               transdata(27*NTEMP + 1:28*NTEMP) = confUSX(:)
               transdata(28*NTEMP + 1:29*NTEMP) = conf2USX(:)
               transdata(29*NTEMP + 1:30*NTEMP) = confstry2(:)
               transdata(30*NTEMP + 1:31*NTEMP) = confstry4(:)
               transdata(31*NTEMP + 1:32*NTEMP) = confUSY(:)
               transdata(32*NTEMP + 1:33*NTEMP) = conf2USY(:)
               transdata(33*NTEMP + 1:34*NTEMP) = cfg_absstr_diag(:)
               transdata(34*NTEMP + 1:35*NTEMP) = cfg_stripediff(:)
               transdata(35*NTEMP + 1:36*NTEMP) = cfg_besten(:)
               transdata(36*NTEMP + 1:37*NTEMP) = cfg_finalen(:)
               transdata(37*NTEMP + 1:38*NTEMP) = confnbar(:)
               transdata(38*NTEMP + 1:39*NTEMP) = conf2nbar(:)
               transdata(39*NTEMP + 1:40*NTEMP) = confnbar2(:)
               transdata(40*NTEMP + 1:41*NTEMP) = confnbar4(:)
               transdata(41*NTEMP + 1:42*NTEMP) = confUBN(:)
               transdata(42*NTEMP + 1:43*NTEMP) = conf2UBN(:)
               transdata(43*NTEMP + 1:44*NTEMP) = cfg_nbar(:)
               transdata(44*NTEMP + 1:45*NTEMP) = cfg_nbar2(:)
               transdata(45*NTEMP + 1:46*NTEMP) = cfg_nbar4(:)
               transdata(46*NTEMP + 1:47*NTEMP) = cfg_nbind(:)

               do itemp = 1, NTEMP
                  transdata(hist_offset + 1 + (itemp - 1)*BINS:hist_offset + itemp*BINS) = confhist(itemp, :)
               end do

               do itemp = 1, NTEMP
                  transdata(flips_offset + 1 + (itemp - 1)*L3:flips_offset + itemp*L3) = &
                  real(cfg_numflips(itemp, :), r8b)
               end do

               transdata(phi_offset + 1:phi_offset + L3) = cfg_phi(:)
               transdata(iconf_offset + 1) = real(iconf, r8b)
               transdata(pt_offset + 1:pt_offset + (NREPLICA - 1)) = real(cfg_ptswap_attempts(:), r8b)
               transdata(pt_offset + NREPLICA:pt_offset + 2*(NREPLICA - 1)) = real(cfg_ptswap_accepts(:), r8b)
               transdata(pt_offset + 2*(NREPLICA - 1) + 1) = real(cfg_ptvisitboth, r8b)
            end if

            call MPI_SEND(transdata, size(transdata), MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
         else
            if (active_config) then
               strxylist(iconf, 1, :) = cfg_strx(:)
               strxylist(iconf, 2, :) = cfg_stry(:)
               absstrlist(iconf, :) = cfg_absstr_diag(:)
               stripedifflist(iconf, :) = cfg_stripediff(:)
               bestenlist(iconf, :) = cfg_besten(:)
               finalenlist(iconf, :) = cfg_finalen(:)
               nbarcfglist(iconf, :) = cfg_nbar(:)
               nbar2cfglist(iconf, :) = cfg_nbar2(:)
               nbar4cfglist(iconf, :) = cfg_nbar4(:)
               nbindcfglist(iconf, :) = cfg_nbind(:)

               summag(:) = summag(:) + confmag(:)
               summag2(:) = summag2(:) + confmag2(:)
               summag4(:) = summag4(:) + confmag4(:)
               sumsusc(:) = sumsusc(:) + confsusc(:)
               sum2susc(:) = sum2susc(:) + conf2susc(:)
               sumUL(:) = sumUL(:) + confUL(:)
               sum2UL(:) = sum2UL(:) + conf2UL(:)
               sum2mag(:) = sum2mag(:) + conf2mag(:)
               sumlogmag(:) = sumlogmag(:) + conflogmag(:)
               sumen(:) = sumen(:) + confen(:)
               sum2en(:) = sum2en(:) + conf2en(:)
               sumsph(:) = sumsph(:) + confsph(:)
               sum2sph(:) = sum2sph(:) + conf2sph(:)
               sumabsstr(:) = sumabsstr(:) + confabsstr(:)
               sumstrx(:) = sumstrx(:) + confstrx(:)
               sumstry(:) = sumstry(:) + confstry(:)
               sumstr2(:) = sumstr2(:) + confstr2(:)
               sumstr4(:) = sumstr4(:) + confstr4(:)
               sumsuscstr(:) = sumsuscstr(:) + confsuscstr(:)
               sum2suscstr(:) = sum2suscstr(:) + conf2suscstr(:)
               sumUS(:) = sumUS(:) + confUS(:)
               sum2US(:) = sum2US(:) + conf2US(:)
               sumnem(:) = sumnem(:) + confnem(:)
               sum2nem(:) = sum2nem(:) + conf2nem(:)
               sum2str(:) = sum2str(:) + conf2str(:)
               sumstrx2(:) = sumstrx2(:) + confstrx2(:)
               sumstrx4(:) = sumstrx4(:) + confstrx4(:)
               sumUSX(:) = sumUSX(:) + confUSX(:)
               sum2USX(:) = sum2USX(:) + conf2USX(:)
               sumstry2(:) = sumstry2(:) + confstry2(:)
               sumstry4(:) = sumstry4(:) + confstry4(:)
               sumUSY(:) = sumUSY(:) + confUSY(:)
               sum2USY(:) = sum2USY(:) + conf2USY(:)
               sumnbar(:) = sumnbar(:) + confnbar(:)
               sum2nbar(:) = sum2nbar(:) + conf2nbar(:)
               sumnbar2(:) = sumnbar2(:) + confnbar2(:)
               sumnbar4(:) = sumnbar4(:) + confnbar4(:)
               sumUBN(:) = sumUBN(:) + confUBN(:)
               sum2UBN(:) = sum2UBN(:) + conf2UBN(:)
               globalhist(:, :) = globalhist(:, :) + confhist(:, :)
               do itemp = 1, NTEMP
                  totalnumflips(iconf, itemp, :) = cfg_numflips(itemp, :)
               end do
               jointp(iconf, :) = cfg_phi(:)
               sum_ptswap_attempts(:) = sum_ptswap_attempts(:) + cfg_ptswap_attempts(:)
               sum_ptswap_accepts(:) = sum_ptswap_accepts(:) + cfg_ptswap_accepts(:)
               sum_ptvisitboth = sum_ptvisitboth + cfg_ptvisitboth
            end if

            do id = 1, numprocs - 1
               call MPI_RECV(transdata, size(transdata), MPI_DOUBLE_PRECISION, id, 1, MPI_COMM_WORLD, status, ierr)
               cfg_id = nint(transdata(iconf_offset + 1))
               if (cfg_id < 1 .or. cfg_id > NCONF) cycle

               strxylist(cfg_id, 1, :) = transdata(14*NTEMP + 1:15*NTEMP)
               strxylist(cfg_id, 2, :) = transdata(15*NTEMP + 1:16*NTEMP)
               absstrlist(cfg_id, :) = transdata(33*NTEMP + 1:34*NTEMP)
               stripedifflist(cfg_id, :) = transdata(34*NTEMP + 1:35*NTEMP)
               bestenlist(cfg_id, :) = transdata(35*NTEMP + 1:36*NTEMP)
               finalenlist(cfg_id, :) = transdata(36*NTEMP + 1:37*NTEMP)
               nbarcfglist(cfg_id, :) = transdata(43*NTEMP + 1:44*NTEMP)
               nbar2cfglist(cfg_id, :) = transdata(44*NTEMP + 1:45*NTEMP)
               nbar4cfglist(cfg_id, :) = transdata(45*NTEMP + 1:46*NTEMP)
               nbindcfglist(cfg_id, :) = transdata(46*NTEMP + 1:47*NTEMP)

               summag(:) = summag(:) + transdata(1:NTEMP)
               summag2(:) = summag2(:) + transdata(NTEMP + 1:2*NTEMP)
               summag4(:) = summag4(:) + transdata(2*NTEMP + 1:3*NTEMP)
               sumsusc(:) = sumsusc(:) + transdata(3*NTEMP + 1:4*NTEMP)
               sum2susc(:) = sum2susc(:) + transdata(4*NTEMP + 1:5*NTEMP)
               sumUL(:) = sumUL(:) + transdata(5*NTEMP + 1:6*NTEMP)
               sum2UL(:) = sum2UL(:) + transdata(6*NTEMP + 1:7*NTEMP)
               sum2mag(:) = sum2mag(:) + transdata(7*NTEMP + 1:8*NTEMP)
               sumlogmag(:) = sumlogmag(:) + transdata(8*NTEMP + 1:9*NTEMP)
               sumen(:) = sumen(:) + transdata(9*NTEMP + 1:10*NTEMP)
               sum2en(:) = sum2en(:) + transdata(10*NTEMP + 1:11*NTEMP)
               sumsph(:) = sumsph(:) + transdata(11*NTEMP + 1:12*NTEMP)
               sum2sph(:) = sum2sph(:) + transdata(12*NTEMP + 1:13*NTEMP)
               sumabsstr(:) = sumabsstr(:) + transdata(13*NTEMP + 1:14*NTEMP)
               sumstrx(:) = sumstrx(:) + transdata(14*NTEMP + 1:15*NTEMP)
               sumstry(:) = sumstry(:) + transdata(15*NTEMP + 1:16*NTEMP)
               sumstr2(:) = sumstr2(:) + transdata(16*NTEMP + 1:17*NTEMP)
               sumstr4(:) = sumstr4(:) + transdata(17*NTEMP + 1:18*NTEMP)
               sumsuscstr(:) = sumsuscstr(:) + transdata(18*NTEMP + 1:19*NTEMP)
               sum2suscstr(:) = sum2suscstr(:) + transdata(19*NTEMP + 1:20*NTEMP)
               sumUS(:) = sumUS(:) + transdata(20*NTEMP + 1:21*NTEMP)
               sum2US(:) = sum2US(:) + transdata(21*NTEMP + 1:22*NTEMP)
               sumnem(:) = sumnem(:) + transdata(22*NTEMP + 1:23*NTEMP)
               sum2nem(:) = sum2nem(:) + transdata(23*NTEMP + 1:24*NTEMP)
               sum2str(:) = sum2str(:) + transdata(24*NTEMP + 1:25*NTEMP)
               sumstrx2(:) = sumstrx2(:) + transdata(25*NTEMP + 1:26*NTEMP)
               sumstrx4(:) = sumstrx4(:) + transdata(26*NTEMP + 1:27*NTEMP)
               sumUSX(:) = sumUSX(:) + transdata(27*NTEMP + 1:28*NTEMP)
               sum2USX(:) = sum2USX(:) + transdata(28*NTEMP + 1:29*NTEMP)
               sumstry2(:) = sumstry2(:) + transdata(29*NTEMP + 1:30*NTEMP)
               sumstry4(:) = sumstry4(:) + transdata(30*NTEMP + 1:31*NTEMP)
               sumUSY(:) = sumUSY(:) + transdata(31*NTEMP + 1:32*NTEMP)
               sum2USY(:) = sum2USY(:) + transdata(32*NTEMP + 1:33*NTEMP)
               sumnbar(:) = sumnbar(:) + transdata(37*NTEMP + 1:38*NTEMP)
               sum2nbar(:) = sum2nbar(:) + transdata(38*NTEMP + 1:39*NTEMP)
               sumnbar2(:) = sumnbar2(:) + transdata(39*NTEMP + 1:40*NTEMP)
               sumnbar4(:) = sumnbar4(:) + transdata(40*NTEMP + 1:41*NTEMP)
               sumUBN(:) = sumUBN(:) + transdata(41*NTEMP + 1:42*NTEMP)
               sum2UBN(:) = sum2UBN(:) + transdata(42*NTEMP + 1:43*NTEMP)
               do itemp = 1, NTEMP
                  totalnumflips(cfg_id, itemp, :) = int(&
     &               transdata(flips_offset + 1 + (itemp - 1)*L3:&
     &               flips_offset + itemp*L3)&
     &            )
               end do
               jointp(cfg_id, :) = transdata(phi_offset + 1:phi_offset + L3)

               do itemp = 1, NTEMP
                  globalhist(itemp, :) = globalhist(itemp, :) + transdata(hist_offset + (itemp - 1)*BINS + 1:hist_offset + itemp*BINS)
               end do
               sum_ptswap_attempts(:) = sum_ptswap_attempts(:) + int(transdata(pt_offset + 1:pt_offset + (NREPLICA - 1)))
               sum_ptswap_accepts(:) = sum_ptswap_accepts(:) + int(transdata(pt_offset + NREPLICA:pt_offset + 2*(NREPLICA - 1)))
               sum_ptvisitboth = sum_ptvisitboth + nint(transdata(pt_offset + 2*(NREPLICA - 1) + 1))

            end do
         end if
#else
         strxylist(iconf, 1, :) = cfg_strx(:)
         strxylist(iconf, 2, :) = cfg_stry(:)
         absstrlist(iconf, :) = cfg_absstr_diag(:)
         stripedifflist(iconf, :) = cfg_stripediff(:)
         bestenlist(iconf, :) = cfg_besten(:)
         finalenlist(iconf, :) = cfg_finalen(:)
         nbarcfglist(iconf, :) = cfg_nbar(:)
         nbar2cfglist(iconf, :) = cfg_nbar2(:)
         nbar4cfglist(iconf, :) = cfg_nbar4(:)
         nbindcfglist(iconf, :) = cfg_nbind(:)

         summag(:) = summag(:) + confmag(:)
         summag2(:) = summag2(:) + confmag2(:)
         summag4(:) = summag4(:) + confmag4(:)
         sumsusc(:) = sumsusc(:) + confsusc(:)
         sum2susc(:) = sum2susc(:) + conf2susc(:)
         sumUL(:) = sumUL(:) + confUL(:)
         sum2UL(:) = sum2UL(:) + conf2UL(:)
         sum2mag(:) = sum2mag(:) + conf2mag(:)
         sumlogmag(:) = sumlogmag(:) + conflogmag(:)
         sumen(:) = sumen(:) + confen(:)
         sum2en(:) = sum2en(:) + conf2en(:)
         sumsph(:) = sumsph(:) + confsph(:)
         sum2sph(:) = sum2sph(:) + conf2sph(:)
         sumabsstr(:) = sumabsstr(:) + confabsstr(:)
         sumstr2(:) = sumstr2(:) + confstr2(:)
         sumstr4(:) = sumstr4(:) + confstr4(:)
         sumsuscstr(:) = sumsuscstr(:) + confsuscstr(:)
         sum2suscstr(:) = sum2suscstr(:) + conf2suscstr(:)
         sumUS(:) = sumUS(:) + confUS(:)
         sum2US(:) = sum2US(:) + conf2US(:)
         sumnem(:) = sumnem(:) + confnem(:)
         sum2nem(:) = sum2nem(:) + conf2nem(:)
         sum2str(:) = sum2str(:) + conf2str(:)
         sumstrx(:) = sumstrx(:) + confstrx(:)
         sumstry(:) = sumstry(:) + confstry(:)
         sumstrx2(:) = sumstrx2(:) + confstrx2(:)
         sumstrx4(:) = sumstrx4(:) + confstrx4(:)
         sumUSX(:) = sumUSX(:) + confUSX(:)
         sum2USX(:) = sum2USX(:) + conf2USX(:)
         sumstry2(:) = sumstry2(:) + confstry2(:)
         sumstry4(:) = sumstry4(:) + confstry4(:)
         sumUSY(:) = sumUSY(:) + confUSY(:)
         sum2USY(:) = sum2USY(:) + conf2USY(:)
         sumnbar(:) = sumnbar(:) + confnbar(:)
         sum2nbar(:) = sum2nbar(:) + conf2nbar(:)
         sumnbar2(:) = sumnbar2(:) + confnbar2(:)
         sumnbar4(:) = sumnbar4(:) + confnbar4(:)
         sumUBN(:) = sumUBN(:) + confUBN(:)
         sum2UBN(:) = sum2UBN(:) + conf2UBN(:)
         do itemp = 1, NTEMP
            totalnumflips(iconf, itemp, :) = cfg_numflips(itemp, :)
         end do
         jointp(iconf, :) = cfg_phi(:)
         sum_ptswap_attempts(:) = sum_ptswap_attempts(:) + cfg_ptswap_attempts(:)
         sum_ptswap_accepts(:) = sum_ptswap_accepts(:) + cfg_ptswap_accepts(:)
         sum_ptvisitboth = sum_ptvisitboth + cfg_ptvisitboth

         globalhist(:, :) = globalhist(:, :) + confhist(:, :)

#endif
         !print *, 'data sent'
      end subroutine Collect_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Write_files(finconf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer(i4b), intent(in) :: finconf
         integer(i4b)   ::  id
         real(r8b)      :: coexprod_sum, coexprod2_sum, coexangle_sum, coexangle2_sum, UBN
!print *, 'DEBUG: WRITESTRIP=', WRITESTRIP

#ifdef PARALLEL
         if (myid == 0) then
#endif
            open (7, file=cfenfile, status='replace')
            rewind (7)
            write (7, *) 'program J1J2_6          '
            write (7, *) 'linear system size      ', L
            write (7, *) 'J1                      ', J1
            write (7, *) 'deltaJ1                 ', deltaJ1
            write (7, *) 'J2                      ', J2
            write (7, *) 'J3                      ', J3
            write (7, *) 'phi, xi                     ', phi, corr_length
            write (7, *) 'impurity conc.          ', PP
            write (7, *) 'impurity mode           ', IMPMODE
            write (7, *) 'start                   ', START
            write (7, *) 'stripe axis             ', STRAXIS
            write (7, *) 'sampler mode            ', SAMPLERMODE
            write (7, *) 'update order            ', UPDATEORDER
            write (7, *) 'temperature protocol    ', TEMPPROTOCOL
            write (7, *) 'PT start mode           ', PTSTARTMODE
            write (7, *) 'PT ladder mode          ', PTLADDERMODE
            write (7, *) 'PT replicas             ', NREPLICA
            write (7, *) 'PT swap every           ', PTSWAPEVER
            write (7, *) 'PT measure index        ', PTMEASINDEX
            write (7, *) 'disorder configurations ', NCONF
            write (7, *) 'equilibration steps     ', NEQ
            write (7, *) 'measurement steps       ', NMESS
            write (7, *) 'RNG seed                ', IRINIT
            if (PTLADDERMODE == 'SPLIT') then
               write (7, *) 'T ladder midpoint       ', TLADDERMID
               write (7, *) 'T ladder low-T fraction ', TLADDERFRAC
            end if
            write (7, *) 'energy definition       ', 'exchange + phi'
            write (7, *) '-------------------------'
            write (7, *) 'disorder configurations processed ', finconf
            write (7, *) '   T      [<energy>]      std.en/sqrt(iconf)      [<spec.heat>]      std.spec.heat/sqrt(iconf)   '
            do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               write (7, '(1x,25(es24.16,1x))') t, sumen(itemp)/finconf, &
                  (sqrt(sum2en(itemp)/finconf - (sumen(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  sumsph(itemp)/finconf, &
                  (sqrt(sum2sph(itemp)/finconf - (sumsph(itemp)/finconf)**2))/sqrt(1.D0*finconf)
            end do
            close (7)

            open (7, file=cfmagfile, status='replace')
            rewind (7)
            write (7, *) 'program J1J2_6          '
            write (7, *) 'spatial system size     ', L
            write (7, *) 'J1                      ', J1
            write (7, *) 'deltaJ1                 ', deltaJ1
            write (7, *) 'J2                      ', J2
            write (7, *) 'J3                      ', J3
            write (7, *) 'phi, xi                     ', phi, corr_length
            write (7, *) 'impurity conc.          ', PP
            write (7, *) 'impurity mode           ', IMPMODE
            write (7, *) 'start                   ', START
            write (7, *) 'stripe axis             ', STRAXIS
            write (7, *) 'sampler mode            ', SAMPLERMODE
            write (7, *) 'update order            ', UPDATEORDER
            write (7, *) 'temperature protocol    ', TEMPPROTOCOL
            write (7, *) 'PT start mode           ', PTSTARTMODE
            write (7, *) 'PT ladder mode          ', PTLADDERMODE
            write (7, *) 'PT replicas             ', NREPLICA
            write (7, *) 'PT swap every           ', PTSWAPEVER
            write (7, *) 'PT measure index        ', PTMEASINDEX
            write (7, *) 'disorder configurations ', NCONF
            write (7, *) 'equilibration steps     ', NEQ
            write (7, *) 'measurement steps       ', NMESS
            write (7, *) 'RNG seed                ', IRINIT
            if (PTLADDERMODE == 'SPLIT') then
               write (7, *) 'T ladder midpoint       ', TLADDERMID
               write (7, *) 'T ladder low-T fraction ', TLADDERFRAC
            end if
            write (7, *) '------------------------'
            write (7, *) 'disorder configurations processed ', finconf
            write (7, *) '   T      [<mag>]      [<mag>^2]      std.mag/sqrt(iconf)      [<susc>]', &
               ' .     std.susc/sqrt(iconf)      [<Binder>]      std.bin/sqrt(iconf) ', &
               '[log <mag>]  Global.Binder'
            do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               write (7, '(1x,25(e12.6,1x))') t, summag(itemp)/finconf, sum2mag(itemp)/finconf, &
                  (sqrt(sum2mag(itemp)/finconf - (summag(itemp)/finconf)**2))/sqrt(1.D0*finconf), sumsusc(itemp)/finconf, &
                  (sqrt(sum2susc(itemp)/finconf - (sumsusc(itemp)/finconf)**2))/sqrt(1.D0*finconf), sumUL(itemp)/finconf, &
                  (sqrt(sum2UL(itemp)/finconf - (sumUL(itemp)/finconf)**2))/sqrt(1.D0*finconf), sumlogmag(itemp)/finconf, &
                  1 - summag4(itemp)/finconf/(3*(summag2(itemp)/finconf)**2)
            end do
            close (7)

            open (7, file=cfnmfile, status='replace')
            rewind (7)
            write (7, *) 'program J1J2_6          '
            write (7, *) 'spatial system size     ', L
            write (7, *) 'J1                      ', J1
            write (7, *) 'deltaJ1                 ', deltaJ1
            write (7, *) 'J2                      ', J2
            write (7, *) 'J3                      ', J3
            write (7, *) 'phi, xi                     ', phi, corr_length
            write (7, *) 'impurity conc.          ', PP
            write (7, *) 'impurity mode           ', IMPMODE
            write (7, *) 'start                   ', START
            write (7, *) 'stripe axis             ', STRAXIS
            write (7, *) 'sampler mode            ', SAMPLERMODE
            write (7, *) 'update order            ', UPDATEORDER
            write (7, *) 'temperature protocol    ', TEMPPROTOCOL
            write (7, *) 'PT start mode           ', PTSTARTMODE
            write (7, *) 'PT ladder mode          ', PTLADDERMODE
            write (7, *) 'PT replicas             ', NREPLICA
            write (7, *) 'PT swap every           ', PTSWAPEVER
            write (7, *) 'PT measure index        ', PTMEASINDEX
            write (7, *) 'disorder configurations ', NCONF
            write (7, *) 'equilibration steps     ', NEQ
            write (7, *) 'measurement steps       ', NMESS
            write (7, *) 'RNG seed                ', IRINIT
            if (PTLADDERMODE == 'SPLIT') then
               write (7, *) 'T ladder midpoint       ', TLADDERMID
               write (7, *) 'T ladder low-T fraction ', TLADDERFRAC
            end if
            write (7, *) '------------------------'
            write (7, *) 'disorder configurations processed ', finconf
            write (7, *) '   T      [<n>]      std.n/sqrt(iconf)      [<n^2>]      [<n^4>]      Global.Nematic.Binder      Nematic.Binder      Nematic.Binder.Error'
            do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               if (sumnbar2(itemp)/finconf > 0.D0) then
                  UBN = 1.D0 - (sumnbar4(itemp)/finconf)/(3.D0*(sumnbar2(itemp)/finconf)**2)
               else
                  UBN = 0.D0
               end if
               write (7, '(1x,25(e12.6,1x))') T, &
                  sumnbar(itemp)/finconf, &
                  (sqrt(sum2nbar(itemp)/finconf - (sumnbar(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  sumnbar2(itemp)/finconf, &
                  sumnbar4(itemp)/finconf, &
                  UBN, &
                  sumUBN(itemp)/finconf, &
                  (sqrt(sum2UBN(itemp)/finconf - (sumUBN(itemp)/finconf)**2))/sqrt(1.D0*finconf)
            end do
            close (7)

            open (7, file=cfstrfile, status='replace')
            rewind (7)
            write (7, *) 'program J1J2_6          '
            write (7, *) 'spatial system size     ', L
            write (7, *) 'J1                      ', J1
            write (7, *) 'deltaJ1                 ', deltaJ1
            write (7, *) 'J2                      ', J2
            write (7, *) 'J3                      ', J3
            write (7, *) 'phi, xi                     ', phi, corr_length
            write (7, *) 'impurity conc.          ', PP
            write (7, *) 'impurity mode           ', IMPMODE
            write (7, *) 'start                   ', START
            write (7, *) 'stripe axis             ', STRAXIS
            write (7, *) 'sampler mode            ', SAMPLERMODE
            write (7, *) 'update order            ', UPDATEORDER
            write (7, *) 'temperature protocol    ', TEMPPROTOCOL
            write (7, *) 'PT start mode           ', PTSTARTMODE
            write (7, *) 'PT ladder mode          ', PTLADDERMODE
            write (7, *) 'PT replicas             ', NREPLICA
            write (7, *) 'PT swap every           ', PTSWAPEVER
            write (7, *) 'PT measure index        ', PTMEASINDEX
            write (7, *) 'disorder configurations ', NCONF
            write (7, *) 'equilibration steps     ', NEQ
            write (7, *) 'measurement steps       ', NMESS
            write (7, *) 'RNG seed                ', IRINIT
            if (PTLADDERMODE == 'SPLIT') then
               write (7, *) 'T ladder midpoint       ', TLADDERMID
               write (7, *) 'T ladder low-T fraction ', TLADDERFRAC
            end if
            write (7, *) '------------------------'
            write (7, *) 'disorder configurations processed ', finconf
            write (7, *) '   T      absstr    str2      str4      nem      Global.Strip.Binder      Strip.susc ', &
               ' Local.Strip.Binder     Local.Strip.Binder.Error      Strip.OP.error      stri.susc.error', &
               ' strx   strx2   strx4   Local.Strx.Binder    Local.Strx.Binder.Error ', &
               ' stry   stry2   stry4   Local.Stry.Binder    Local.Stry.Binder.Error', &
               ' absstrx_absstry  absstrx_absstry.Error  coex_angle  coex_angle.Error'
            do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               coexprod_sum = 0.D0
               coexprod2_sum = 0.D0
               coexangle_sum = 0.D0
               coexangle2_sum = 0.D0
               do id = 1, finconf
                  coexprod_sum = coexprod_sum + strxylist(id, 1, itemp)*strxylist(id, 2, itemp)
                  coexangle_sum = coexangle_sum + Stripe_coexistence_angle(strxylist(id, 1, itemp), strxylist(id, 2, itemp))
                  coexprod2_sum = coexprod2_sum + (strxylist(id, 1, itemp)*strxylist(id, 2, itemp))**2
                  coexangle2_sum = coexangle2_sum + Stripe_coexistence_angle(strxylist(id, 1, itemp), strxylist(id, 2, itemp))**2
               end do
               write (7, '(1x,25(e12.6,1x))') T, &
                  sumabsstr(itemp)/finconf, &
                  sumstr2(itemp)/finconf, &
                  sumstr4(itemp)/finconf, &
                  sumnem(itemp)/finconf, &
                  (2 - sumstr4(itemp)/finconf/(sumstr2(itemp)/finconf)**2), &
                  sumsuscstr(itemp)/finconf, &
                  sumUS(itemp)/finconf, &
                  (sqrt(sum2US(itemp)/finconf - (sumUS(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  (sqrt(sum2str(itemp)/finconf - (sumabsstr(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  (sqrt(sum2suscstr(itemp)/finconf - (sumsuscstr(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  sumstrx(itemp)/finconf, &
                  sumstrx2(itemp)/finconf, &
                  sumstrx4(itemp)/finconf, &
                  sumUSX(itemp)/finconf, &
                  (sqrt(sum2USX(itemp)/finconf - (sumUSX(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  sumstry(itemp)/finconf, &
                  sumstry2(itemp)/finconf, &
                  sumstry4(itemp)/finconf, &
                  sumUSY(itemp)/finconf, &
                  (sqrt(sum2USY(itemp)/finconf - (sumUSY(itemp)/finconf)**2))/sqrt(1.D0*finconf), &
                  coexprod_sum/finconf, &
                  (sqrt(coexprod2_sum/finconf - (coexprod_sum/finconf)**2))/sqrt(1.D0*finconf), &
                  coexangle_sum/finconf, &
                  (sqrt(coexangle2_sum/finconf - (coexangle_sum/finconf)**2))/sqrt(1.D0*finconf)
            end do
            close (7)

            if (WRITEHIST .eqv. .true.) then
               open (7, file=cfhsfile, status='replace')
               rewind (7)
               write (7, *) 'program J1J2_6          '
               write (7, *) 'linear system size      ', L
               write (7, *) 'J1                      ', J1
               write (7, *) 'deltaJ1                 ', deltaJ1
               write (7, *) 'J2                      ', J2
               write (7, *) 'J3                      ', J3
               write (7, *) 'phi                     ', phi
               write (7, *) 'impurity conc.          ', PP
               write (7, *) 'impurity mode           ', IMPMODE
               write (7, *) 'start                   ', START
               write (7, *) 'stripe axis             ', STRAXIS
               write (7, *) 'sampler mode            ', SAMPLERMODE
               write (7, *) 'update order            ', UPDATEORDER
               write (7, *) 'temperature protocol    ', TEMPPROTOCOL
               write (7, *) 'PT start mode           ', PTSTARTMODE
               write (7, *) 'PT ladder mode          ', PTLADDERMODE
               write (7, *) 'PT replicas             ', NREPLICA
               write (7, *) 'PT swap every           ', PTSWAPEVER
               write (7, *) 'PT measure index        ', PTMEASINDEX
               write (7, *) 'disorder configurations ', NCONF
               write (7, *) 'equilibration steps     ', NEQ
               write (7, *) 'measurement steps       ', NMESS
               write (7, *) 'RNG seed                ', IRINIT
               if (PTLADDERMODE == 'SPLIT') then
                  write (7, *) 'T ladder midpoint       ', TLADDERMID
                  write (7, *) 'T ladder low-T fraction ', TLADDERFRAC
               end if
               write (7, *) 'energy definition       ', 'exchange + phi'
               write (7, *) '-------------------------'
               write (7, *) 'disorder configurations processed ', finconf
               write (7, *) 'T  BIN   [<energy>]  freq'

               ! Write histogram bins for this temperature
               do itemp = 1, NTEMP
               T = Temperature_for_index(itemp)
               do ibin = 1, BINS - 1
                  write (7, '(1x, 1(e12.6, 1x), I4, 1x, 2(e12.6, 1x))') T, ibin, &
                     hist_en_min_t(itemp) + hist_bin_size_t(itemp)*(ibin - 1), &
                     globalhist(itemp, ibin)/(finconf)
               end do
               end do

            end if

            close (7)
            if (WRITESTRIP) then

               open (7, file=cfxyfile, status='replace')
               rewind (7)
               write (7, *) 'iconf, T, strx, stry, absstr, D, best_en, final_en, absstrx_absstry, coex_angle'

               do id = 1, finconf
                  do itemp = 1, NTEMP
                     T = Temperature_for_index(itemp)
                     write (7, '(I4,1x,9(E12.6,1x))') id, T, strxylist(id, 1, itemp), strxylist(id, 2, itemp), &
                        absstrlist(id, itemp), stripedifflist(id, itemp), bestenlist(id, itemp), finalenlist(id, itemp), &
                        strxylist(id, 1, itemp)*strxylist(id, 2, itemp), &
                        Stripe_coexistence_angle(strxylist(id, 1, itemp), strxylist(id, 2, itemp))
                  end do
               end do
            end if
            close (7)

            open (7, file=cfndfile, status='replace')
            rewind (7)
            write (7, *) 'iconf, T, nbar, nbar2, nbar4, nbinder'
            do id = 1, finconf
               do itemp = 1, NTEMP
                  T = Temperature_for_index(itemp)
                  write (7, '(I4,1x,5(E12.6,1x))') id, T, nbarcfglist(id, itemp), nbar2cfglist(id, itemp), &
                     nbar4cfglist(id, itemp), nbindcfglist(id, itemp)
               end do
            end do
            close (7)

            if (WRITENUMFLIPS) then

               open (7, file=cfnumfile, status='replace')
               rewind (7)
               write (7, *) 'iconf, T, isite, numflips'

               do id = 1, finconf
                  do itemp = 1, NTEMP
                     T = Temperature_for_index(itemp)
                     do isite = 0, L3 - 1
                        write (7, '(I4,1x,E12.6,1x,I8,1x,I8)') id, T, isite, totalnumflips(id, itemp, isite)
                     end do
                  end do
               end do
            end if
            close (7)

            if (SAMPLERMODE == 'PT') then
               open (7, file=cfptfile, status='replace')
               rewind (7)
               write (7, *) 'pair, T_low, T_high, attempts, accepts, acceptance'
               do itemp = 1, NREPLICA - 1
                  if (sum_ptswap_attempts(itemp) > 0) then
                     write (7, '(I4,1x,2(E12.6,1x),2(I12,1x),E12.6)') itemp, temp_ladder(itemp), temp_ladder(itemp + 1), &
                        sum_ptswap_attempts(itemp), sum_ptswap_accepts(itemp), &
                        real(sum_ptswap_accepts(itemp), r8b)/real(sum_ptswap_attempts(itemp), r8b)
                  else
                     write (7, '(I4,1x,2(E12.6,1x),2(I12,1x),E12.6)') itemp, temp_ladder(itemp), temp_ladder(itemp + 1), &
                        0, 0, 0.D0
                  end if
               end do
               write (7, *) 'replicas_visiting_both_extremes ', sum_ptvisitboth, &
                  real(sum_ptvisitboth, r8b)/real(max(1, finconf*NREPLICA), r8b)
               close (7)
            end if

            if (WRITEPHITABLE) then

               open (7, file=cfptable, status='replace')
               rewind (7)
               write (7, *) 'iconf, isite, p'

               do id = 1, finconf

                  do isite = 0, L3 - 1
                     write (7, '(I4,1x,I8,1x,E12.6)') id, isite, jointp(id, isite)

                  end do
               end do
            end if
            close (7)

#ifdef PARALLEL
         end if           ! of if (myid==0)
#endif

      end subroutine Write_files

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
!     print *,c1
         if (c1 .ne. 536870911) then
            print *, 'Nonstandard integer representation. Stoped.'
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
