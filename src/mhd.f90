program mhd
  !-- Resout du/dt + df/dx = 0 pour la MHD 1D (test de Brio-Wu)
  !-- Compilation : gfortran -fdefault-real-8 mhd_solvers.f90 mhd.f90 -o mhd
  !-- Execution   : .\mhd
  !-- Trace       : python mhd.py

  use mhd_solvers
  implicit none

  integer, parameter :: nx = 800   !-- Nombre de cellules
  integer, parameter :: ngc = 2    !-- Cellules fantomes
  real, parameter :: cfl = 0.9     !-- Coefficient CFL
  real, parameter :: xmax = 1.0    !-- Domaine [0,xmax]
  real, parameter :: tmax = 0.1    !-- Temps final

  real :: x(nx)
  real :: v(nvar, 1-ngc:nx+ngc)    !-- Primitives (avec cellules fantomes)
  real :: u(nvar, nx), u0(nvar, nx)
  real :: t, dt, dx
  integer :: i, k

  !-- Grille
  dx = xmax / real(nx)
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do

  !-- Conditions initiales de Brio-Wu : (rho, vx, vy, vz, By, Bz, p)
  do i = 1, nx
    if (x(i) <= 0.5) then
      v(:,i) = [ 1.0,   0.0, 0.0, 0.0,  1.0, 0.0, 1.0 ]
    else
      v(:,i) = [ 0.125, 0.0, 0.0, 0.0, -1.0, 0.0, 0.1 ]
    end if
    call prim_to_cons(v(:,i), u(:,i))
  end do

  t = 0.0

  !-- Evolution : Runge-Kutta d'ordre 2 (Heun)
  do while (t < tmax)
    u0 = u
    call evolve(0.0, 1.0, .true.)     ! pas d'Euler, calcule dt
    call evolve(0.5, 0.5, .false.)    ! correcteur de Heun, meme dt
    t = t + dt
    print *, 't =', t
  end do

  !-- Ecriture : x puis les 7 primitives
  open(10, file = 'mhd.txt')
  do i = 1, nx
    write(10, '(8(es16.6,1x))') x(i), v(:,i)
  end do
  close(10)

contains

!-------------------------------------

subroutine evolve(c1, c2, set_dt)
  real, intent(in) :: c1, c2
  logical, intent(in) :: set_dt
  real :: vl(nvar), vr(nvar), f(nvar, 0:nx), amax, a
  integer :: i, k

  !-- Conditions aux limites (extrapolation d'ordre zero)
  do k = 1, nvar
    v(k, -1)   = v(k, 1)
    v(k, 0)    = v(k, 1)
    v(k, nx+1) = v(k, nx)
    v(k, nx+2) = v(k, nx)
  end do

  amax = 0.0

  !-- Flux a toutes les interfaces i+1/2
  do i = 0, nx
    do k = 1, nvar
      vl(k) = reconstruction( v(k,i-1), v(k,i),   v(k,i+1) )
      vr(k) = reconstruction( v(k,i+2), v(k,i+1), v(k,i)   )
    end do
    call llf(vl, vr, f(:,i), a)
    amax = max(amax, a)
  end do

  !-- Pas de temps (CFL), fixe au premier sous-pas seulement
  if (set_dt) then
    dt = cfl * dx / amax
    if (t + dt > tmax) dt = tmax - t    ! tombe exactement sur tmax
  end if

  !-- Mise a jour conservative puis retour aux primitives
  do i = 1, nx
    u(:,i) = c1 * u0(:,i) + c2 * ( u(:,i) - (dt/dx) * (f(:,i) - f(:,i-1)) )
    call cons_to_prim(u(:,i), v(:,i))
  end do

end subroutine evolve

end program mhd