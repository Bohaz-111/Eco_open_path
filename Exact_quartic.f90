program Exact_quartic
   use iso_fortran_env, only: real64
   implicit none

   ! Compile with LAPACK and BLAS:
   ! gfortran -O2 Exact_quartic.f90 -llapack -lblas -o Exact_quartic
   ! ./Exact_quartic

   integer,  parameter :: dp = real64
   integer,  parameter :: nbasis = 10           ! Basis states n = 0, ..., nbasis-1
   integer,  parameter :: ngrid = 401
   real(dp), parameter :: pi = 3.14159265358979323846_dp
   real(dp), parameter :: temp = 0.125_dp
   real(dp), parameter :: beta = 1.0_dp/temp
   real(dp), parameter :: pmax = 4.0_dp
   complex(dp), parameter :: iu = (0.0_dp, 1.0_dp) ! Imaginary unit

   ! V(x) = x**2/2 + x**4/4; m = hbar = k_B = 1.
   ! Increase nbasis to check basis convergence, pmax to check momentum tails.

   real(dp) :: h(0:nbasis-1,0:nbasis-1), p2(0:nbasis-1,0:nbasis-1)
   real(dp) :: energy(0:nbasis-1), pop(0:nbasis-1)
   real(dp) :: work(3*nbasis), pgrid(ngrid), np(ngrid)
   real(dp) :: integral, mean_p2
   complex(dp) :: chi(0:nbasis-1), psi(0:nbasis-1)
   integer :: i, n, info

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

   call build_hamiltonian(h, p2)
   call dsyev('V', 'U', nbasis, h, nbasis, energy, work, size(work), info)
   if (info /= 0) error stop 'Hamiltonian diagonalisation failed'
   ! energy(j) = E_j and h(n,j) = <n|psi_j>, with n,j = 0, ..., nbasis-1.
   ! The ground state is energy(0), with coefficients h(:,0).

   pop = exp(-beta*(energy-energy(0)))
   pop = pop/sum(pop)
   mean_p2 = sum(pop*sum(h*matmul(p2,h), dim=1)) ! Calculate <p^2> for later checking. 

   ! Main calculation of n(p) starts here. 
   do i = 1, ngrid
      pgrid(i) = pmax*real(i-1,dp)/real(ngrid-1,dp)

      ! Momentum-space basis: chi(n) = <p|n>.
      chi(0) = pi**(-0.25_dp)*exp(-0.5_dp*pgrid(i)**2) ! chi(0) = <p|0>, see P5 - Eq.(7)
      chi(1) = -iu*sqrt(2.0_dp)*pgrid(i)*chi(0) ! chi(1) = <p|1>, see the text below P5 - Eq.(7)

      ! Implements the recurrence relation P5 - Eq.(6).
      do n = 1, nbasis-2
         chi(n+1) = sqrt(real(n,dp)/real(n+1,dp))*chi(n-1) &
                  - iu*pgrid(i)*sqrt(2.0_dp/real(n+1,dp))*chi(n)
      end do

      psi = matmul(chi, h) ! psi(j) = <p|psi_j> = sum_n <p|n><n|psi_j>
      np(i) = sum(pop*abs(psi)**2) ! n(p) = sum_j w_j |psi(j)|^2
   end do

   call write_two_col('exact_momentum_10basis.dat', pgrid, np)
   integral = 2.0_dp*pmax/real(ngrid-1,dp)*(sum(np)-0.5_dp*(np(1)+np(ngrid))) ! Check normalisation. 

   print *, 'Ground-state energy:', energy(0)
   print *, 'Thermal energy:', sum(pop*energy)
   print *, '<p**2> from the operator:', mean_p2
   print *, 'n(0):', np(1)
   print *, '2 * integral from 0 to pmax of n(p) dp (target 1):', integral

contains

   subroutine build_hamiltonian(h, p2)
      real(dp), intent(out) :: h(0:,0:), p2(0:,0:)
      integer, parameter :: npad = nbasis+4 ! Allow for 4 sequential raising/lowering ladder operators.
      real(dp) :: x(0:npad-1,0:npad-1), x2(0:npad-1,0:npad-1), x4(0:npad-1,0:npad-1)
      integer :: n

      ! Unit-frequency oscillator basis. Pad BEFORE multiplying operators so
      ! paths through states outside the retained basis contribute to x**4.
      x = 0.0_dp
      do n = 0, npad-2
         x(n,n+1) = sqrt(real(n+1,dp)/2.0_dp)
         x(n+1,n) = x(n,n+1)
      end do
      x2 = matmul(x,x)
      x4 = matmul(x2,x2)

      ! Since H_0 = a^dagger a + 1/2 = (p^2+x^2)/2, p^2 = 2a^dagger a + 1 - x^2.
      ! Matrix elements: <n|p^2|m> = -<n|x^2|m> + (2n+1) delta_{nm}.
      p2 = -x2(0:nbasis-1,0:nbasis-1)
      do n = 0, nbasis-1
         p2(n,n) = p2(n,n)+real(2*n+1,dp)
      end do

      h = 0.5_dp*p2 + 0.5_dp*x2(0:nbasis-1,0:nbasis-1) + 0.25_dp*x4(0:nbasis-1,0:nbasis-1)
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
