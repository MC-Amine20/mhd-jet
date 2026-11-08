program jet2d
  !-- Jet hydrodynamique 2D (tranche) se propageant selon x
  !-- Compilation : gfortran -fdefault-real-8 hydro2d_solvers.f90 jet2d.f90 -o jet2d

  use hydro2d_solvers
  implicit none

  integer, parameter :: nx = 240, ny = 160, ngc = 2
  real, parameter :: cfl = 0.4
  real, parameter :: xmax = 12.0, ymax = 8.0
  real, parameter :: tmax = 2.4
  integer, parameter :: nout = 24         !-- nombre de snapshots

  !-- Milieu ambiant et jet
  real, parameter :: rho_a = 1.0, p_a = 1.0
  real, parameter :: eta = 0.1, mach = 5.0, rj = 1.0
  real, parameter :: pfloor = 1.0e-6, rhofloor = 1.0e-6

  real :: x(nx), y(ny)
  real :: v(nvar, 1-ngc:nx+ngc, 1-ngc:ny+ngc)
  real :: u(nvar, nx, ny), u0(nvar, nx, ny)
  real :: vjet(nvar)
  real :: t, dt, dx, dy, yc, rho_j, p_j, vj, dtout, tnext
  integer :: i, j, n, nstep

  dx = xmax / real(nx)
  dy = ymax / real(ny)
  yc = 0.5 * ymax
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do
  do j = 1, ny
    y(j) = (real(j) - 0.5) * dy
  end do

  !-- Etat du jet : (rho, vx, vy, p)
  rho_j = eta * rho_a
  p_j   = p_a
  vj    = mach * sqrt(gamma * p_j / rho_j)
  vjet  = [ rho_j, vj, 0.0, p_j ]
  print *, 'Vitesse du jet vj =', vj
  print *, 'Vitesse theorique de la tete (1D) ~', vj*sqrt(eta)/(1.0+sqrt(eta))

  !-- Conditions initiales : milieu ambiant au repos
  do j = 1, ny
    do i = 1, nx
      v(:,i,j) = [ rho_a, 0.0, 0.0, p_a ]
      call prim_to_cons(v(:,i,j), u(:,i,j))
    end do
  end do

  t = 0.0
  dtout = tmax / real(nout)
  tnext = dtout
  n = 0
  nstep = 0

  do while (n < nout)
    u0 = u
    call evolve(0.0, 1.0, .true.)
    call evolve(0.5, 0.5, .false.)
    t = t + dt
    nstep = nstep + 1
    if (mod(nstep, 200) == 0) print *, 'pas', nstep, ' t =', t
    if (t >= tnext - 1.0e-10) then
      n = n + 1
      call write_snapshot(n)
      tnext = tnext + dtout
    end if
  end do

contains

!-------------------------------------

subroutine write_snapshot(num)
  integer, intent(in) :: num
  character(len=20) :: fname
  integer :: i, j

  write(fname, '(a,i3.3,a)') 'jet_', num, '.txt'
  open(10, file = fname)
  do i = 1, nx
    do j = 1, ny
      write(10, '(6(es14.5,1x))') x(i), y(j), v(1,i,j), v(2,i,j), v(3,i,j), v(4,i,j)
    end do
  end do
  close(10)
end subroutine write_snapshot

!-------------------------------------

subroutine evolve(c1, c2, set_dt)
  real, intent(in) :: c1, c2
  logical, intent(in) :: set_dt
  real :: vl(nvar), vr(nvar), a, amaxx, amaxy
  real :: fx(nvar, 0:nx, ny), fy(nvar, nx, 0:ny)
  integer :: i, j, k

  !-- Bord gauche (x=0) : injection du jet ou mur reflechissant
  do j = 1, ny
    if (abs(y(j) - yc) <= rj) then
      v(:, 0,  j) = vjet
      v(:, -1, j) = vjet
    else
      v(:, 0,  j) = v(:, 1, j);  v(2, 0,  j) = -v(2, 1, j)
      v(:, -1, j) = v(:, 2, j);  v(2, -1, j) = -v(2, 2, j)
    end if
  end do

  !-- Bord droit : sortie libre
  do j = 1, ny
    v(:, nx+1, j) = v(:, nx, j)
    v(:, nx+2, j) = v(:, nx, j)
  end do

  !-- Bords haut et bas : sortie libre
  do i = -1, nx+2
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

  !-- Pas de temps (on tombe exactement sur chaque snapshot)
  if (set_dt) then
    dt = cfl / ( amaxx/dx + amaxy/dy )
    if (t + dt > tnext) dt = tnext - t
  end if

  !-- Mise a jour conservative + planchers de securite
  do j = 1, ny
    do i = 1, nx
      u(:,i,j) = c1 * u0(:,i,j) + c2 * ( u(:,i,j) &
               - (dt/dx) * (fx(:,i,j) - fx(:,i-1,j)) &
               - (dt/dy) * (fy(:,i,j) - fy(:,i,j-1)) )
      call cons_to_prim(u(:,i,j), v(:,i,j))
      if (v(1,i,j) < rhofloor .or. v(4,i,j) < pfloor) then
        v(1,i,j) = max(v(1,i,j), rhofloor)
        v(4,i,j) = max(v(4,i,j), pfloor)
        call prim_to_cons(v(:,i,j), u(:,i,j))
      end if
    end do
  end do

end subroutine evolve

end program jet2d