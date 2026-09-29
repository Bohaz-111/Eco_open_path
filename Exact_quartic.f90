program Exact_quartic
   use iso_fortran_env, only: real64
   implicit none

   ! Compile with LAPACK and BLAS:
   ! gfortran -O2 Exact_quartic.f90 -llapack -lblas -o Exact_quartic
   ! ./Exact_quartic

   integer,  parameter :: dp = real64
   integer,  parameter :: nbasis = 60           ! Harmonic oscillator basis states, 5 is enough, we used 60 here is more than needed.
   integer,  parameter :: ngrid = 401
   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: temp = 0.125_dp
   real(dp), parameter :: beta = 1.0_dp/temp
   real(dp), parameter :: pmax = 4.0_dp
   real(dp), parameter :: c2 = 1.0_dp, c4 = 1.0_dp

   ! V(x) = c2*x**2/2 + c4*x**4/4; m = hbar = k_B = 1.
   ! Increase nbasis to check basis convergence, pmax to check momentum tails.

   real(dp) :: h(nbasis,nbasis), p2(nbasis,nbasis), energy(nbasis), pop(nbasis)
   real(dp) :: work(3*nbasis), phi(nbasis), pgrid(ngrid), np(ngrid)
   real(dp) :: integral, mean_p2
   complex(dp) :: phase(nbasis), psi(nbasis)
   integer :: i, k, info

   interface
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*), work(*)
         integer, intent(out) :: info
      end subroutine dsyev
   end interface

   if (nbasis < 2 .or. ngrid < 2 .or. temp <= 0.0_dp .or. pmax <= 0.0_dp) error stop 'Invalid parameters'
   if (c4 < 0.0_dp .or. (c4 <= 0.0_dp .and. c2 <= 0.0_dp)) error stop 'Potential must be confining'

   call build_hamiltonian(h, p2)
   call dsyev('V', 'U', nbasis, h, nbasis, energy, work, size(work), info)
   if (info /= 0) error stop 'Hamiltonian diagonalisation failed'
   ! Now h(:,k) contains the expansion coefficients of energy eigenstate k.

   pop = exp(-beta*(energy-energy(1)))
   pop = pop/sum(pop)
   mean_p2 = sum(pop*sum(h*matmul(p2,h), dim=1))

   ! Fourier transform of oscillator basis state n is (-i)**n times itself.
   do k = 1, nbasis
      phase(k) = cmplx(0.0_dp, -1.0_dp, kind=dp)**(k-1)
   end do

   do i = 1, ngrid
      pgrid(i) = pmax*real(i-1,dp)/real(ngrid-1,dp)

      ! Normalised Hermite functions, with oscillator quantum number n = k-1.
      phi(1) = pi**(-0.25_dp)*exp(-0.5_dp*pgrid(i)**2)
      phi(2) = sqrt(2.0_dp)*pgrid(i)*phi(1)
      do k = 2, nbasis-1
         phi(k+1) = sqrt(2.0_dp/real(k,dp))*pgrid(i)*phi(k) &
                  - sqrt(real(k-1,dp)/real(k,dp))*phi(k-1)
      end do

      psi = matmul(phi*phase, h)
      np(i) = sum(pop*abs(psi)**2)
   end do

   call write_two_col('exact_momentum_60basis.dat', pgrid, np)
   integral = 2.0_dp*pmax/real(ngrid-1,dp)*(sum(np)-0.5_dp*(np(1)+np(ngrid)))

   print *, 'Ground-state energy:', energy(1)
   print *, 'Thermal energy:', sum(pop*energy)
   print *, '<p**2> from the operator:', mean_p2
   print *, 'n(0):', np(1)
   print *, '2 * integral from 0 to pmax of n(p) dp (target 1):', integral
   ! No finite-window renormalisation: omitted tails account for missing area.

contains

   subroutine build_hamiltonian(h, p2)
      real(dp), intent(out) :: h(:,:), p2(:,:)
      integer, parameter :: n = nbasis+4
      real(dp) :: x(n,n), x2(n,n), x4(n,n)
      integer :: k

      ! Unit-frequency oscillator basis. Pad BEFORE multiplying operators so
      ! paths through states outside the retained basis contribute to x**4.
      x = 0.0_dp
      do k = 1, n-1
         x(k,k+1) = sqrt(real(k,dp)/2.0_dp)
         x(k+1,k) = x(k,k+1)
      end do
      x2 = matmul(x,x)
      x4 = matmul(x2,x2)

      p2 = -x2(1:nbasis,1:nbasis)
      do k = 1, nbasis
         p2(k,k) = p2(k,k)+real(2*k-1,dp)
      end do

      h = 0.5_dp*p2 + 0.5_dp*c2*x2(1:nbasis,1:nbasis) + 0.25_dp*c4*x4(1:nbasis,1:nbasis)
   end subroutine build_hamiltonian

   subroutine write_two_col(filename, x, y)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: x(:), y(:)
      integer :: i

      ! Columns: p >= 0, n(p). Full-line density; no factor 2 or 4*pi*p**2.
      open(unit=20, file=filename, status='replace', action='write')
      do i = 1, size(x)
         write(20,'(2(ES24.15E3,1X))') x(i), y(i)
      end do
      close(20)
   end subroutine write_two_col

end program Exact_quartic
