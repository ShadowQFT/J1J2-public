module corr_rand_3d_mod
   implicit none
   private
   public :: make_corr_rand_3d
contains
!=============================================================
! 3D correlated random field via Fourier filtering
!
! Steps:
!   1) Build real-space correlation C(dx,dy,dz) = exp(-r/xi)
!      with periodic distances.
!   2) FFT C -> Ck
!   3) Generate white Gaussian noise X
!   4) FFT X -> Xk
!   5) Xk <- Xk * sqrt(|Ck|)
!   6) Inverse FFT -> correlated Gaussian field
!   7) Optional map Gaussian -> Uniform in [-MAX_VALUE, +MAX_VALUE]
!
! Output: i  j  k  value  to ising_corr_rand_3d.dat
!
! Compile:
!   gfortran -O2 -Wall ising_corr_rand_3d.f90 fft2d/fftsg3d.f fft2d/fftsg.f -o ising3d
! Run:
!   ./ising3d
!=============================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates complementary error function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION erfcc(x)
   implicit none
   integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)   ! 8-byte reals
   real(r8b)   :: erfcc, x
   real(r8b)   :: t, z
   z = abs(x)
   t = 1.D0/(1.D0 + 0.5D0*z)
   erfcc = t*exp(-z*z - 1.26551223D0 + t*(1.00002368D0 + t*(.37409196D0 + t*&
  &(.09678418D0 + t*(-.18628806D0 + t*(.27886807D0 + t*(-1.13520398D0 + t*&
  &(1.48851587D0 + t*(-.82215223D0 + t*.17087277D0)))))))))
   if (x .lt. 0.D0) erfcc = 2.D0 - erfcc
   return
END FUNCTION erfcc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gaussian random number generator gkiss05
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION gkiss05()
   implicit none
   integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)
   integer, parameter      :: i4b = SELECTED_INT_KIND(8)
   real(r8b)             :: gkiss05
   real(r8b), external   :: rkiss05
   real(r8b)             :: v1, v2, s, fac
   integer(i4b)          :: iset
   real(r8b)             :: gset
   common/gausscom/gset, iset
   if (iset .ne. 1) then
      do
         v1 = 2.D0*rkiss05() - 1.D0
         v2 = 2.D0*rkiss05() - 1.D0
         s = v1*v1 + v2*v2
         if ((s < 1.D0) .and. (s > 0.D0)) exit
      end do
      fac = sqrt(-2.D0*log(s)/s)
      gset = v1*fac
      iset = 1
      gkiss05 = v2*fac
   else
      iset = 0
      gkiss05 = gset
   end if
   return
END FUNCTION gkiss05

SUBROUTINE gkissinit(iinit)
   implicit none
   integer, parameter     :: r8b = SELECTED_REAL_KIND(P=14, R=99)
   integer, parameter     :: i4b = SELECTED_INT_KIND(8)
   integer(i4b)          :: iinit, iset
   real(r8b)             :: gset
   common/gausscom/gset, iset
   iset = 0
   call kissinit(iinit)
end subroutine gkissinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KISS05 random number generator (uniform [0,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION rkiss05()
   implicit none
   integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)
   integer, parameter      :: i4b = SELECTED_INT_KIND(8)
   real(r8b), parameter    :: am = 4.656612873077392578d-10
   real(r8b)             :: rkiss05
   integer(i4b)          :: kiss
   integer(i4b)          :: x, y, z, w
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
   integer, parameter      :: r8b = SELECTED_REAL_KIND(P=14, R=99)
   integer, parameter     :: i4b = SELECTED_INT_KIND(8)
   integer(i4b) idum, ia, im, iq, ir, iinit
   integer(i4b) k, x, y, z, w, c1, c2, c3, c4
   real(r8b) rkiss05, rdum
   parameter(ia=16807, im=2147483647, iq=127773, ir=2836)
   common/kisscom/x, y, z, w
   c1 = -8
   c1 = ishftc(c1, -3)
   if (c1 .ne. 536870911) then
      print *, 'Nonstandard integer representation. Stoped.'
      stop
   end if
   idum = iinit
   idum = abs(1099087573*idum)
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

!=============================================================
! NOTE FOR COPYING INTO YOUR MAIN SIM FILE
!
! If you want ZERO extra source files:
!   1) Copy/paste `make_corr_rand_3d` (above) into your simulation.
!   2) Also copy/paste the helper routines below (erfcc, gkiss05,
!      gkissinit, rkiss05, kissinit) into the same simulation file
!      OR replace the gkiss05/gkissinit calls with your own RNG.
!
! You will STILL need to link the Ooura FFT sources that provide:
!   - rdft3d
!   - rdft3dsort
!=============================================================

!=============================================================
! 3D correlated random field via Fourier filtering (Ooura FFT)
!
! Copy/paste `make_corr_rand_3d` into your main simulation file.
!
! Requirements in your build (still required):
!   - Link Ooura FFT routines that provide `rdft3d` and `rdft3dsort`
!   - Provide RNG helpers `gkissinit`, `gkiss05` and mapping helper `erfcc`
!     (they are included later in this file and can also be copy/pasted).
!
! Output:
!   x(0:n1-1,0:n2-1,0:n3-1) correlated Gaussian field with approx
!   C(r) ~ exp(-r/xi) using periodic distances.
!=============================================================

subroutine make_corr_rand_3d(x, n1, n2, n3, xi, seed, do_uniform, max_value)
   implicit none
   integer, parameter :: r8b = selected_real_kind(14, 99)
   integer, parameter :: i4b = selected_int_kind(8)

   integer, intent(in) :: n1, n2, n3
   real(r8b), intent(in) :: xi
   integer(i4b), intent(in) :: seed
   logical, intent(in) :: do_uniform
   real(r8b), intent(in) :: max_value

   real(r8b), intent(out) :: x(0:n1-1, 0:n2-1, 0:n3-1)

   ! FFTSG requires power-of-two lengths; we pad internally.
   integer :: n1fft, n2fft, n3fft
   ! Ooura RDFT3D packing requires first dimension >= n1fft+2
   integer :: n1max, n2max

   ! Work sizes
   integer :: nwork, ip_len, w_len, t_len

   ! Main working arrays (include padding/packing in first dim)
   real(r8b), allocatable :: xw(:,:,:)
   real(r8b), allocatable :: c(:,:,:)

   ! FFT work arrays
   real(r8b), allocatable :: t(:)
   real(r8b), allocatable :: w(:)
   integer,  allocatable :: ip(:)

   ! Loop indices
   integer :: i, j, k
   integer :: k1, k2, k3

   ! External FFT routines (from fftsg3d.f)
   external :: rdft3d
   external :: rdft3dsort

   ! Local scaling
   real(r8b) :: scale

   !--------------------------------------------------------
   ! Derived sizes / allocate
   !--------------------------------------------------------
   n1fft = nextpow2(n1)
   n2fft = nextpow2(n2)
   n3fft = nextpow2(n3)

   n1max = n1fft + 2
   n2max = n2fft

   nwork  = max(n1fft/2, max(n2fft, n3fft))
   ip_len = 2 + int(sqrt(dble(nwork))) + 10
   w_len  = max(n1fft/4, max(n2fft/2, n3fft/2)) + n1fft/4 + 1024
   t_len  = max(8*n2fft, 8*n3fft)

   allocate(xw(0:n1max-1, 0:n2max-1, 0:n3fft-1))
   allocate(c(0:n1max-1, 0:n2max-1, 0:n3fft-1))
   allocate(t(0:t_len-1))
   allocate(w(0:w_len-1))
   allocate(ip(0:ip_len-1))

   !--------------------------------------------------------
   ! Seed RNG
   !--------------------------------------------------------
   call gkissinit(seed)

   !--------------------------------------------------------
   ! Build real-space correlation C(dx,dy,dz)
   !--------------------------------------------------------
   c(:,:,:) = 0.0D0
   do k = 0, n3fft - 1
      do j = 0, n2fft - 1
         do i = 0, n1fft - 1
            c(i, j, k) = corr3d_exp(i, j, k, n1fft, n2fft, n3fft, xi)
         end do
      end do
   end do

   !--------------------------------------------------------
   ! FFT correlation -> Ck
   !--------------------------------------------------------
   ip(0) = 0
   call rdft3d(n1max, n2max, n1fft, n2fft, n3fft, 1, c, t, ip, w)
   call rdft3dsort(n1max, n2max, n1fft, n2fft, n3fft, 1, c)

   !--------------------------------------------------------
   ! Generate white Gaussian noise X
   !--------------------------------------------------------
   xw(:,:,:) = 0.0D0
   do k = 0, n3fft - 1
      do j = 0, n2fft - 1
         do i = 0, n1fft - 1
            xw(i, j, k) = gkiss05()
         end do
      end do
   end do

   !--------------------------------------------------------
   ! FFT noise -> Xk
   !--------------------------------------------------------
   call rdft3d(n1max, n2max, n1fft, n2fft, n3fft, 1, xw, t, ip, w)
   call rdft3dsort(n1max, n2max, n1fft, n2fft, n3fft, 1, xw)

   !--------------------------------------------------------
   ! Fourier-space filtering: Xk <- Xk * sqrt(|Ck|)
   !
   ! After sort:
   !   xw(2*k1,  k2, k3) = Re
   !   xw(2*k1+1,k2, k3) = Im
   ! for 0 <= k1 <= n1/2
   !--------------------------------------------------------
   do k3 = 0, n3fft - 1
      do k2 = 0, n2fft - 1
         do k1 = 0, n1fft/2
            xw(2*k1,     k2, k3) = xw(2*k1,     k2, k3) * sqrt(abs(c(2*k1, k2, k3)))
            xw(2*k1 + 1, k2, k3) = xw(2*k1 + 1, k2, k3) * sqrt(abs(c(2*k1, k2, k3)))
         end do
      end do
   end do

   !--------------------------------------------------------
   ! Inverse FFT back to real space
   !--------------------------------------------------------
   call rdft3dsort(n1max, n2max, n1fft, n2fft, n3fft, -1, xw)
   call rdft3d(n1max, n2max, n1fft, n2fft, n3fft, -1, xw, t, ip, w)

   ! Ooura scaling (matches your original code)
   scale = 2.0D0 / dble(n1fft) / dble(n2fft) / dble(n3fft)

   ! Copy out the physical region and scale
   do k = 0, n3 - 1
      do j = 0, n2 - 1
         do i = 0, n1 - 1
            x(i, j, k) = xw(i, j, k) * scale
         end do
      end do
   end do

   !--------------------------------------------------------
   ! Optional: Gaussian -> Uniform in [-max_value, +max_value]
   !--------------------------------------------------------
   if (do_uniform) then
      do k = 0, n3 - 1
         do j = 0, n2 - 1
            do i = 0, n1 - 1
               x(i, j, k) = (2.0D0*(1.0D0 - 0.5D0*erfcc(x(i, j, k)/sqrt(2.0D0))) - 1.0D0) * max_value
            end do
         end do
      end do
   end if

   deallocate(xw, c, t, w, ip)

contains

   !-----------------------------------------------------------
   ! Small helper: next power of two >= n
   !-----------------------------------------------------------
   integer function nextpow2(n) result(p)
      implicit none
      integer, intent(in) :: n
      p = 1
      do while (p < n)
         p = p * 2
      end do
   end function nextpow2

   !-----------------------------------------------------------
   ! Periodic-distance exponential correlation
   !-----------------------------------------------------------
   real(r8b) function corr3d_exp(ix, iy, iz, nx, ny, nz, xi) result(corr)
      implicit none
      integer, intent(in) :: ix, iy, iz, nx, ny, nz
      real(r8b), intent(in) :: xi
      integer :: dx, dy, dz
      real(r8b) :: r

      dx = min(ix, nx - ix)
      dy = min(iy, ny - iy)
      dz = min(iz, nz - iz)

      r = sqrt(dble(dx*dx + dy*dy + dz*dz))

      if (xi <= 0.0D0) then
         if (r == 0.0D0) then
            corr = 1.0D0
         else
            corr = 0.0D0
         end if
      else
         corr = exp(-r/xi)
      end if
   end function corr3d_exp

end subroutine make_corr_rand_3d
end module corr_rand_3d_mod
