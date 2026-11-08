program hydro2d
  !-- Euler 2D, test de Sod le long de x (idir=1) ou de y (idir=2)
  !-- Compilation : gfortran -fdefault-real-8 hydro2d_solvers.f90 hydro2d.f90 -o hydro2d

  use hydro2d_solvers
  implicit none

  integer, parameter :: idir = 1       !-- 1 : Sod selon x, 2 : Sod selon y
  integer, parameter :: nx = 400       !-- mettre 4 si idir = 2
  integer, parameter :: ny = 4         !-- mettre 400 si idir = 2
  integer, parameter :: ngc = 2
  real, parameter :: cfl = 0.4
  real, parameter :: xmax = 1.0, ymax = 1.0
  real, parameter :: tmax = 0.2

  real :: x(nx), y(ny)
  real :: v(nvar, 1-ngc:nx+ngc, 1-ngc:ny+ngc)
  real :: u(nvar, nx, ny), u0(nvar, nx, ny)
  real :: t, dt, dx, dy, s
  integer :: i, j

  dx = xmax / real(nx)
  dy = ymax / real(ny)
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do
  do j = 1, ny
    y(j) = (real(j) - 0.5) * dy
  end do

  !-- Conditions initiales : (rho, vx, vy, p)
  do j = 1, ny
    do i = 1, nx
      if (idir == 1) then
        s = x(i)
      else
        s = y(j)
      end if
      if (s <= 0.5) then
        v(:,i,j) = [ 1.0,   0.0, 0.0, 1.0 ]
      else
        v(:,i,j) = [ 0.125, 0.0, 0.0, 0.1 ]
      end if
      call prim_to_cons(v(:,i,j), u(:,i,j))
    end do
  end do

  t = 0.0
  do while (t < tmax)
    u0 = u
    call evolve(0.0, 1.0, .true.)
    call evolve(0.5, 0.5, .false.)
    t = t + dt
    print *, 't =', t
  end do

  !-- Sortie : profil le long de la direction du choc (x, rho, v normale, p)
  open(10, file = 'hydro2d.txt')
  if (idir == 1) then
    do i = 1, nx
      write(10, '(4(es16.6,1x))') x(i), v(1,i,1), v(2,i,1), v(4,i,1)
    end do
  else
    do j = 1, ny
      write(10, '(4(es16.6,1x))') y(j), v(1,1,j), v(3,1,j), v(4,1,j)
    end do
  end if
  close(10)

contains

!-------------------------------------

subroutine evolve(c1, c2, set_dt)
  real, intent(in) :: c1, c2
  logical, intent(in) :: set_dt
  real :: vl(nvar), vr(nvar), a, amaxx, amaxy
  real :: fx(nvar, 0:nx, ny), fy(nvar, nx, 0:ny)
  integer :: i, j, k

  !-- Conditions aux limites (sortie libre : copie de la derniere cellule)
  do j = 1, ny
    v(:, 0,    j) = v(:, 1,  j)
    v(:, -1,   j) = v(:, 1,  j)
    v(:, nx+1, j) = v(:, nx, j)
    v(:, nx+2, j) = v(:, nx, j)
  end do
  do i = 1, nx
    v(:, i, 0   ) = v(:, i, 1 )
    v(:, i, -1  ) = v(:, i, 1 )
    v(:, i, ny+1) = v(:, i, ny)
    v(:, i, ny+2) = v(:, i, ny)
  end do

  amaxx = 0.0
  amaxy = 0.0

  !-- Flux selon x aux interfaces (i+1/2, j)
  do j = 1, ny
    do i = 0, nx
      do k = 1, nvar
        vl(k) = reconstruction( v(k,i-1,j), v(k,i,  j), v(k,i+1,j) )
        vr(k) = reconstruction( v(k,i+2,j), v(k,i+1,j), v(k,i,  j) )
      end do
      call llf(vl, vr, 1, fx(:,i,j), a)
      amaxx = max(amaxx, a)
    end do
  end do

  !-- Flux selon y aux interfaces (i, j+1/2)
  do j = 0, ny
    do i = 1, nx
      do k = 1, nvar
        vl(k) = reconstruction( v(k,i,j-1), v(k,i,j  ), v(k,i,j+1) )
        vr(k) = reconstruction( v(k,i,j+2), v(k,i,j+1), v(k,i,j  ) )
      end do
      call llf(vl, vr, 2, fy(:,i,j), a)
      amaxy = max(amaxy, a)
    end do
  end do

  !-- Pas de temps (CFL en 2D)
  if (set_dt) then
    dt = cfl / ( amaxx/dx + amaxy/dy )
    if (t + dt > tmax) dt = tmax - t
  end if

  !-- Mise a jour conservative
  do j = 1, ny
    do i = 1, nx
      u(:,i,j) = c1 * u0(:,i,j) + c2 * ( u(:,i,j) &
               - (dt/dx) * (fx(:,i,j) - fx(:,i-1,j)) &
               - (dt/dy) * (fy(:,i,j) - fy(:,i,j-1)) )
      call cons_to_prim(u(:,i,j), v(:,i,j))
    end do
  end do

end subroutine evolve

end program hydro2d