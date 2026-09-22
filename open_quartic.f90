program open_quartic
   use iso_fortran_env, only: real64
   implicit none

   ! To compile:
   ! gfortran -O3 open_quartic.f90 -o open_quartic
   ! ./open_quartic

   integer,  parameter :: dp    = real64
   integer,  parameter :: P     = 64              ! Number of potential beads
   real(dp), parameter :: pi    = 3.14159265358979323846_dp
   real(dp), parameter :: temp  = 0.125_dp
   real(dp), parameter :: tau   = 5.0_dp          ! Bussi thermostat time
   real(dp), parameter :: beta  = 1.0_dp/temp
   real(dp), parameter :: beta_n = beta/real(P, dp)
   real(dp), parameter :: dt    = 0.1_dp*beta_n
   integer,  parameter :: steps = 1000000
   integer,  parameter :: equil_steps = 50000
   integer,  parameter :: ntraj = 1
   integer,  parameter :: stride = 10
   integer,  parameter :: ngrid = 401
   real(dp), parameter :: pmax  = 4.0_dp

   real(dp) :: q(P), mom(P), f(P)
   real(dp) :: pgrid(ngrid), np(ngrid), cos_sum(ngrid)
   real(dp) :: delta, closure_sum, integral
   integer  :: i, traj

   ! V(x)=x**4/4; m=hbar=k_B=1.
   ! P-1 springs, no closing spring.
   call random_seed()
   do i = 1, ngrid
      pgrid(i) = pmax*real(i-1, dp)/real(ngrid-1, dp)
   end do
   cos_sum = 0.0_dp
   closure_sum = 0.0_dp

   do traj = 1, ntraj
      print *, 'Trajectory', traj, '/', ntraj
      call init_config(q, mom)
      do i = 1, equil_steps
         call step_vv(q, mom, f, dt)
         call thermostat(mom, dt, tau, 1.0_dp/beta_n)
      end do

      do i = 1, steps
         call step_vv(q, mom, f, dt)
         call thermostat(mom, dt, tau, 1.0_dp/beta_n)
         if (mod(i, stride) /= 0) cycle
         delta = q(P) - q(1)
         cos_sum = cos_sum + cos(pgrid*delta)
         closure_sum = closure_sum + exp(-delta**2/(2.0_dp*beta_n))
      end do
   end do

   np = sqrt(beta_n/(2.0_dp*pi))*exp(-0.5_dp*beta_n*pgrid**2) * cos_sum/closure_sum
   ! Columns: p >= 0, n(p). In 1D there is no angular Jacobian.
   call write_two_col('momentum.dat', pgrid, np, ngrid)
   integral = (pmax/real(ngrid-1, dp)) * (sum(np) - 0.5_dp*(np(1) + np(ngrid)))
   ! By symmetry, twice the positive-half integral should approach one.
   print *, '2 * integral from 0 to pmax of n(p) dp (target 1):', 2.0_dp*integral

contains

   subroutine init_config(q, mom)
      real(dp), intent(out) :: q(:), mom(:)
      integer :: j

      q = 0.0_dp
      do j = 1, size(mom)
         mom(j) = sqrt(1.0_dp/beta_n)*randn()
      end do
   end subroutine init_config

   pure subroutine force(q, f)
      real(dp), intent(in)  :: q(:)
      real(dp), intent(out) :: f(:)
      real(dp) :: spring
      integer :: j

      f = -q**3
      do j = 1, size(q)-1
         spring = (q(j+1) - q(j))/beta_n**2
         f(j)   = f(j)   + spring
         f(j+1) = f(j+1) - spring
      end do
   end subroutine force

   pure subroutine step_vv(q, mom, f, dt)
      real(dp), intent(inout) :: q(:), mom(:), f(:)
      real(dp), intent(in)    :: dt
      real(dp) :: halfdt

      halfdt = 0.5_dp*dt
      call force(q, f)
      mom = mom + halfdt*f
      q   = q   + dt*mom
      call force(q, f)
      mom = mom + halfdt*f
   end subroutine step_vv

   subroutine thermostat(mom, dt, tau, T) ! Bussi/CSVR thermostat, 2007 JCP
      real(dp), intent(inout) :: mom(:)
      real(dp), intent(in)    :: dt, tau, T
      integer  :: dof, kshape
      real(dp) :: K, Kbar, c, s, r1, sum_r2, alpha2, alpha, factor

      dof  = size(mom)
      K    = 0.5_dp*sum(mom*mom)
      Kbar = 0.5_dp*real(dof, dp)*T
      c = exp(-dt/tau)
      s = 1.0_dp-c

      r1 = randn()
      ! Chi-square draw with dof-1 degrees of freedom, for even or odd P.
      kshape = (dof-1)/2
      sum_r2 = rand_gamma(kshape)
      if (mod(dof-1, 2) == 1) sum_r2 = sum_r2 + randn()**2

      factor = Kbar/(real(dof, dp)*max(K, tiny(1.0_dp)))
      alpha2 = c + factor*s*(r1*r1 + sum_r2) &
               + 2.0_dp*exp(-0.5_dp*dt/tau)*sqrt(factor*s)*r1
      alpha = sign(sqrt(max(alpha2, 0.0_dp)), sqrt(c) + sqrt(factor*s)*r1)
      mom = alpha*mom
   end subroutine thermostat

   real(dp) function rand_gamma(k)
      integer, intent(in) :: k
      real(dp) :: v(k)

      call random_number(v)
      v = max(v, 1.0e-12_dp)
      rand_gamma = -2.0_dp*sum(log(v))
   end function rand_gamma

   real(dp) function randn()
      real(dp) :: u1, u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      randn = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function randn

   subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(dp),     intent(in) :: x(:), y(:)
      integer,      intent(in) :: n
      integer :: u, i

      open(newunit=u, file=fname, status='replace', action='write')
      do i = 1, n
         write(u, '(ES24.15E3,1X,ES24.15E3)') x(i), y(i)
      end do
      close(u)
      print *, 'Saved to: ', fname
   end subroutine write_two_col

end program open_quartic
