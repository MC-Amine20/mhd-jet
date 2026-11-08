program hydro
  !-- Solves du/dt + df/dx = 0 for 1D Euler equations
  !-- grid x, primitive variables v, conservative variables u (and u0)
  !-- To compile:  gfortran hydro_solvers.f90 hydro.f90 -o hydro
  !-- To run:      ./hydro
  !-- To plot:     python hydro.py

  use hydro_solvers
  implicit none

  integer, parameter :: nx = 400   !-- Number of grid points
  integer, parameter :: ngc = 2    !-- Number of ghost cells
  real, parameter :: cfl = 0.9     !-- CFL coefficient (cfl<1)
  real, parameter :: xmax = 1.0    !-- Spatial domain [0,xmax]
  real, parameter :: tmax = 0.2    !-- Output time

  real :: x(nx)                    !-- Grid
  real :: v(nvar, 1-ngc:nx+ngc)    !-- Primitive variables (with ghost cells)
  real :: u(nvar, nx), u0(nvar, nx)  !-- Conservative variables (u0 also needed)
  real :: t, dt, dx
  integer :: i, k

  !-- Grid definition
  dx = xmax / real(nx)
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do

  !-- Initial conditions on v, then v -> u
  do i = 1, nx
    if (x(i) <= 0.5) then
      v(:,i) = [ 1.0, 0.0, 1.0 ]           ! left state
    else
      v(:,i) = [ 0.125, 0.0, 0.1 ]         ! right state
    end if
    call prim_to_cons(v(:,i), u(:,i))
  end do

  !-- Initialize ghost cells (extrapolation)
  do k = 1, nvar
    v(k, 0)  = v(k, 1)
    v(k, -1) = v(k, 1)
    v(k, nx+1) = v(k, nx)
    v(k, nx+2) = v(k, nx)
  end do

  t = 0.0

  !-- Time evolution using 2nd order Runge-Kutta (Heun)
  do while (t < tmax)
    u0 = u
    call evolve(0.0, 1.0)    ! Euler step
    call evolve(0.5, 0.5)    ! Heun corrector
    t = t + dt
    print *, 't =', t        ! monitor progress
  end do

  !-- Write results to text file
  open(10, file = 'hydro.txt')
  do i = 1, nx
    write(10, '(4(es16.6,1x))') x(i), v(1,i), v(2,i), v(3,i)
  end do
  close(10)

contains

!-------------------------------------

subroutine evolve(c1, c2)
  !-- Update u and v for each RK substep
  real, intent(in) :: c1, c2   !-- RK coefficients
  real :: vl(nvar), vr(nvar), f(nvar, 0:nx), amax, a
  integer :: i, k

  !-- Boundary conditions (zero‑order extrapolation)
  do k = 1, nvar
    v(k, -1) = v(k, 1)
    v(k, 0)  = v(k, 1)
    v(k, nx+1) = v(k, nx)
    v(k, nx+2) = v(k, nx)
  end do

  amax = 0.0

  !-- Compute fluxes at all interfaces i+1/2 (i=0..nx)
  do i = 0, nx
    do k = 1, nvar
      vl(k) = reconstruction( v(k,i-1), v(k,i), v(k,i+1) ) ! left state
      vr(k) = reconstruction( v(k,i+2), v(k,i+1), v(k,i) ) ! right state
    end do
    call llf(vl, vr, f(:,i), a)
    amax = max(amax, a)
  end do

  !-- Time step (CFL condition)
  dt = cfl * dx / amax
  if (dt <= 0.0) then
    print *, 'ERROR: non‑positive dt'
    stop
  end if

  !-- Conservative update and conversion back to primitive
  do i = 1, nx
    u(:,i) = c1 * u0(:,i) + c2 * ( u(:,i) - (dt/dx) * (f(:,i) - f(:,i-1)) )
    call cons_to_prim(u(:,i), v(:,i))
  end do

end subroutine evolve

!-------------------------------------

end program hydro