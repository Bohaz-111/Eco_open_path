program Eco_open_quartic
   use iso_fortran_env, only: real64
   implicit none

   ! To compile:
   ! gfortran -O3 Eco_open_quartic.f90 -o Eco_open_quartic
   ! ./Eco_open_quartic

   integer,  parameter :: dp    = real64

   integer,  parameter :: mode  = 1               ! 0 = Trotter, 1 = Eco
   integer,  parameter :: P     = 32              ! Number of potential beads

   real(dp), parameter :: pi    = 3.14159265358979323846_dp
   real(dp), parameter :: temp  = 0.125_dp
   real(dp), parameter :: tau   = 5.0_dp          ! Bussi thermostat time
   real(dp), parameter :: beta  = 1.0_dp/temp
   real(dp), parameter :: beta_n = beta/real(P, dp)
   real(dp), parameter :: dt    = 0.1_dp*beta_n

   integer,  parameter :: steps = 1000000
   integer,  parameter :: equil_steps = 50000
   integer,  parameter :: ntraj = 10
   integer,  parameter :: stride = 10

   integer,  parameter :: ngrid = 401
   real(dp), parameter :: pmax  = 4.0_dp

   real(dp), parameter :: omega_max = 4.0_dp      ! Eco harmonic fitting window
   integer,  parameter :: nfit = max(100, ceiling(10.0_dp*beta*omega_max))

   real(dp) :: q(P), mom(P), f(P)
   real(dp) :: pgrid(ngrid), np(ngrid), cos_sum(ngrid)
   real(dp) :: delta, closure_sum, integral

   real(dp) :: spring_correction(P,P), sigma2
   real(dp) :: fit_x(nfit), fit_t(nfit), fit_g(nfit), weight(P-1)

   integer  :: i, traj

   ! V(x)=x**2/2+x**4/4; m=hbar=k_B=1.
   ! Trotter: P-1 springs, no closing spring. Eco adds a fitted spring correction.
   ! Original joint fit: equal fractional errors in endpoint variance and open gyration.

   if (mode /= 0 .and. mode /= 1) stop 'mode must be 0 (Trotter) or 1 (Eco)'

   sigma2 = beta_n
   if (mode == 1) call init_eco()
   call random_seed()

   do i = 1, ngrid
      pgrid(i) = pmax*real(i-1, dp)/real(ngrid-1, dp)
   end do

   cos_sum = 0.0_dp
   closure_sum = 0.0_dp

   do traj = 1, ntraj
      print *, 'Trajectory', traj, '/', ntraj
      call init_config(q, mom)

      do i = 1, equil_steps ! Equilibration in NVT.
         if (mode == 1) mom = mom - 0.5_dp*dt*matmul(spring_correction, q)
         call step_vv(q, mom, f, dt)
         if (mode == 1) mom = mom - 0.5_dp*dt*matmul(spring_correction, q)
         call thermostat(mom, dt, tau, 1.0_dp/beta_n)
      end do

      do i = 1, steps       ! Production in NVT. 
         if (mode == 1) mom = mom - 0.5_dp*dt*matmul(spring_correction, q)
         call step_vv(q, mom, f, dt)
         if (mode == 1) mom = mom - 0.5_dp*dt*matmul(spring_correction, q)
         call thermostat(mom, dt, tau, 1.0_dp/beta_n)
         if (mod(i, stride) /= 0) cycle

         delta = q(P) - q(1)
         cos_sum = cos_sum + cos(pgrid*delta)
         closure_sum = closure_sum + exp(-delta**2/(2.0_dp*sigma2))
      end do
   end do

   np = sqrt(sigma2/(2.0_dp*pi))*exp(-0.5_dp*sigma2*pgrid**2) * cos_sum/closure_sum

   ! Columns: p >= 0, n(p). In 1D there is no angular Jacobian.
   if (mode == 0) then
      call write_two_col('Trotter_momentum.dat', pgrid, np, ngrid)
   else
      call write_two_col('EcoP32_momentum.dat', pgrid, np, ngrid)
   end if

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

      f = -q-q**3

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

   subroutine thermostat(mom, dt, tau, T) ! Bussi/CSVR thermostat, 2007 JCP.
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

   subroutine init_eco()
      real(dp) :: y(P-1), basis(P), z, u, rms(2), r(2*nfit)
      ! y(k) = beta*Omega_k
      ! basis(j) is component j of normal mode basis vector.
      ! z is temporary scalar
      ! u is dimensionless kernel variance, sigma^2/beta
      ! rms(1) is endpoint RMS fractional error.
      ! rms(2) is gyration RMS fractional error.
      ! r is two blocks of scaled fractional errors.
      real(dp), save :: jac(2*nfit,P-1) ! Jacobian, sensitivity of error with respect to parameters.
      integer :: j, k

      if (P < 1 .or. temp <= 0.0_dp .or. omega_max <= 0.0_dp) stop 'Invalid Eco parameters'

      z = beta*omega_max/(2.0_dp*real(P,dp))
      u = tanh(z)/(z*real(P,dp))
      sigma2 = beta*u

      ! x = beta*omega; targets are endpoint variance/beta and open gyration/beta.
      do j = 1, nfit
         fit_x(j) = beta*omega_max*(real(j,dp)-0.5_dp)/real(nfit,dp)
         z = fit_x(j)
         fit_t(j) = 2.0_dp*tanh(0.5_dp*z)/z ! Harmonic target end to end distance.
         fit_g(j) = (z/tanh(z)-1.0_dp)/(2.0_dp*z**2) ! Harmonic target radius of gyration.
      end do

      weight = 0.0_dp
      do k = 1, P-1, 2
         weight(k) = 8.0_dp*cos(pi*real(k,dp)/(2.0_dp*P))**2 ! Weights only on odd modes.
      end do

      if (P > 1) call fit_modes(y)
      call fit_residual(log(y), r, jac) ! Calculates fractional error and slope of error due to changing parameters.
      rms = sqrt([sum(r(1:nfit)**2), sum(r(nfit+1:)**2)])

      print *, 'Eco fit RMS fractional errors (endpoint, gyration):', rms
      print *, 'Endpoint Gaussian variance:', sigma2, ' (Trotter:', beta_n, ')'
      if (maxval(rms) > 0.05_dp) print *, 'Warning: coarse harmonic fit; increase P.' ! Warn if RMS fractional error > 5%

      spring_correction = 0.0_dp
      do k = 1, P-1
         do j = 1, P
            basis(j) = sqrt(2.0_dp/real(P,dp))*cos(pi*(real(j,dp)-0.5_dp)*real(k,dp)/P)
         end do

         do j = 1, P
            spring_correction(:,j) = spring_correction(:,j)+(y(k)/beta)**2*basis*basis(j)
         end do
      end do

      ! K_Eco - K_Trotter: correction kicks preserve the original force/VV routines.
      do j = 1, P-1
         spring_correction(j,j) = spring_correction(j,j) - 1.0_dp/beta_n**2
         spring_correction(j+1,j+1) = spring_correction(j+1,j+1) - 1.0_dp/beta_n**2
         spring_correction(j,j+1) = spring_correction(j,j+1) + 1.0_dp/beta_n**2
         spring_correction(j+1,j) = spring_correction(j+1,j) + 1.0_dp/beta_n**2
      end do

      if (P > 1) print *, 'Largest fitted spring frequency:', maxval(y)/beta

      ! Columns: cosine normal-mode index k, Trotter frequency, fitted Eco frequency.
      open(unit=20, file='Eco_modes.dat', status='replace', action='write')
      write(20,'(I6,2(1X,ES24.15E3))') 0, 0.0_dp, 0.0_dp
      do k = 1, P-1
         write(20,'(I6,2(1X,ES24.15E3))') k, &
            2.0_dp*sin(pi*real(k,dp)/(2.0_dp*real(P,dp)))/beta_n, y(k)/beta
      end do
      close(20)
   end subroutine init_eco

   subroutine fit_residual(z, r, jac)
      ! Equal-weight fractional errors; z_k = log(beta*omega_k) ensures positivity in omega_k.
      real(dp), intent(in) :: z(:) ! Takes in adjustable parameters (log internal polymer frequencies)
      real(dp), intent(out) :: r(:), jac(:,:) ! Output fractional error
                                              ! and sensitivity of this error due to changing the parameters Jac(i, k) = ∂r(i)/∂z(k).
      real(dp) :: y2(size(z)), inv(size(z)), deriv(size(z)), scale
      integer :: j

      y2 = exp(2.0_dp*z) ! = (βΩ_k)^2
      scale = sqrt(real(nfit,dp))

      do j = 1, nfit
         inv = 1.0_dp/(fit_x(j)**2+y2) ! = d_{jk}
         deriv = -2.0_dp*y2*inv**2 ! = ∂d_{jk}/∂z_k

         r(j) = ((sum(weight*inv)+sigma2/beta)/fit_t(j)-1.0_dp)/scale ! End to end distance fractional error/sqrt{N}.
         r(nfit+j) = (sum(inv)/fit_g(j)-1.0_dp)/scale ! Radius of gyration fractional error/sqrt{N}. 

         jac(j,:) = weight*deriv/(fit_t(j)*scale) ! = J_{ik} = ∂r_i/∂z_k for end to end distance.
         jac(nfit+j,:) = deriv/(fit_g(j)*scale) ! = J_{ik} = ∂r_i/∂z_k for radius of gyration.
      end do
   end subroutine fit_residual

   subroutine fit_modes(y)
      ! Levenberg-Marquardt, P4
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      real(dp), intent(out) :: y(:)
      real(dp) :: z(size(y)), trial(size(y)), grad(size(y)), projected(size(y))
      real(dp) :: r(2*nfit), rt(2*nfit), jac(2*nfit,size(y)), jt(2*nfit,size(y))
      real(dp) :: a(size(y),size(y)), damping, cost, best_cost, lower, upper, pg
      integer, parameter :: max_iter = 800
      integer :: seed, iter, k
      logical :: ok, converged, best_converged
      logical :: active(size(y))
      character(len=9) :: start_name
      character(len=20) :: status

      lower = log(1.0e-6_dp)
      upper = log(4.0_dp*real(P,dp))
      best_cost = huge(1.0_dp)
      best_converged = .false.

      do seed = 1, 2
         start_name = 'Trotter'
         if (seed == 2) start_name = 'Matsubara'
         do k = 1, size(y)
            z(k) = 2.0_dp*P*sin(pi*real(k,dp)/(2.0_dp*P)) ! Dimensionless Trotter frequency beta*Omega_k.
            if (seed == 2) z(k) = pi*real(k,dp) ! Matsubara frequency start. 
         end do
         z = log(z)

         damping = 1.0e-3_dp ! Lambda in P4.
         converged = .false.
         status = 'iteration limit'
         call fit_residual(z, r, jac) ! Calculate the residue/error and Jacobian of starting frequencies.
         cost = 0.5_dp*sum(r*r) ! = C = Mean square error. 

         ! iter counts attempted steps; also check the final accepted point.
         do iter = 0, max_iter
            pg = huge(1.0_dp)
            if (.not. ieee_is_finite(cost) .or. .not. all(ieee_is_finite(jac))) then
               status = 'non-finite fit'
               exit
            end if
            grad = matmul(transpose(jac),r) ! = g_k = ∂C/∂z_k Implements P4 - Eq.(4)
            if (.not. all(ieee_is_finite(grad))) then
               status = 'non-finite gradient'
               exit
            end if

            ! Identify components whose downhill components are blocked.
            active = (z >= upper-1.0e-10_dp .and. grad < 0.0_dp) .or. &
                     (z <= lower+1.0e-10_dp .and. grad > 0.0_dp)
            projected = grad
            where (active) projected = 0.0_dp
            pg = maxval(abs(projected))

            if (pg < 1.0e-9_dp) then
               converged = .true.
               status = 'converged'
               exit
            end if
            if (iter == max_iter) exit
            if (damping > 1.0e12_dp) then
               status = 'damping limit'
               exit
            end if

            a = matmul(transpose(jac),jac)
            do k = 1, size(y)
               a(k,k) = a(k,k)+damping*max(a(k,k),1.0e-10_dp) ! Implements P4 - Eq.(11), enforce positive definiteness to set up for Cholesky.
            end do                                            ! A is also real and symmetric. 

            ! Freeze active variables in BOTH the matrix and right-hand side.
            ! Clipping a fully coupled step afterwards can turn it uphill.
            do k = 1, size(y)
               if (active(k)) then
                  a(k,:) = 0.0_dp
                  a(:,k) = 0.0_dp
                  a(k,k) = 1.0_dp
               end if
            end do

            call solve_positive(a, -projected, trial, ok)
            ! Solve A δz = -g, to find δz using Cholesky, ok reports whether the solve has suceeded.

            if (ok) then
               trial = min(upper,max(lower,z+trial))
               call fit_residual(trial, rt, jt) ! Calculate the residue and Jacobian of candidate.

               ! Same cost decrease, with less cancellation than subtracting two costs.
               if (all(ieee_is_finite(rt)) .and. all(ieee_is_finite(jt))) then
                  if (0.5_dp*sum((r-rt)*(r+rt)) > 0.0_dp) then
                     z = trial
                     r = rt
                     jac = jt
                     cost = 0.5_dp*sum(r*r)
                     damping = max(1.0e-12_dp,damping/3.0_dp)
                     cycle
                  end if
               end if
            end if

            damping = damping*10.0_dp ! We get here if either Cholesky solve failed or candiate has bigger cost.
                                      ! Tread more carefully. (i.e. bigger damping)
         end do

         write(*,'(A,A,A,I0,A,ES12.4,A,ES12.4,2A)') 'Eco fit ',trim(start_name),': steps=',iter, &
              ', cost=',cost,', projected gradient=',pg,', ',trim(status)

         if (converged .and. cost < best_cost) then
            best_cost = cost
            y = exp(z)
            best_converged = .true.
         end if
      end do

      if (.not. best_converged) error stop 'Eco fit failed: neither start met the gradient tolerance.'
   end subroutine fit_modes

   subroutine solve_positive(a, b, x, ok)
      ! Cholesky solve for the damped least-squares step
      ! Solve A x = b
      ! Factorise A = L L^T, where L is a lower trangular matrix and L^T is its transpose. 
      ! Then L L^T x = L y = b, solve for y first using forward substitution.
      ! Then L^T x = y, finally solve for x using back substitution.
      real(dp), intent(inout) :: a(:,:) ! Modified during factorisation, lower triangle eventually overwritten by L. 
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp) :: s
      integer :: i, j, n

      n = size(b)
      ok = .false.

      do i = 1, n
         do j = 1, i ! The double loop only visits the lower triangle.
            s = a(i,j)-dot_product(a(i,1:j-1),a(j,1:j-1)) ! = L(i, j)*L(j, j)
            ! Not obvious, but a(i, 1:j-1) is actually L(i, 1:j-1) since it a has been overwritten.
            if (i == j) then
               if (s <= 0.0_dp) return ! i.e. The subroutine does not work with a non positive-definte matrix. 
               a(i,j) = sqrt(s) ! Implements P3 - Eq.(6)
            else
               a(i,j) = s/a(j,j) ! Implements P3 - Eq.(5)
            end if
         end do
      end do

      ! Forward substitution
      do i = 1, n
         x(i) = (b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/a(i,i) ! Implements P3 - Eq.(8)
      end do

      ! Back substitution
      do i = n, 1, -1
         x(i) = (x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/a(i,i) ! Implements P3 - Eq.(9)
      end do

      ok = .true.
   end subroutine solve_positive

end program Eco_open_quartic
