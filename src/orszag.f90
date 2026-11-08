program orszag
  !-- Vortex d'Orszag-Tang (MHD 2D, nettoyage GLM de div B), domaine periodique
  !-- Compilation : gfortran -O2 -fdefault-real-8 mhd2d_solvers.f90 orszag.f90 -o orszag

  use mhd2d_solvers
  implicit none

  integer, parameter :: nx = 128, ny = 128, ngc = 2
  real, parameter :: cfl = 0.4
  real, parameter :: xmax = 1.0, ymax = 1.0
  real, parameter :: tmax = 0.5
  real, parameter :: alpha_glm = 0.1       !-- amortissement de psi
  real, parameter :: pi = 3.14159265358979323846

  real :: x(nx), y(ny)
  real :: v(nvar, 1-ngc:nx+ngc, 1-ngc:ny+ngc)
  real :: u(nvar, nx, ny), u0(nvar, nx, ny)
  real :: t, dt, dx, dy, rho0, p0, b0
  integer :: i, j, nstep

  dx = xmax / real(nx)
  dy = ymax / real(ny)
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do
  do j = 1, ny
    y(j) = (real(j) - 0.5) * dy
  end do

  rho0 = 25.0 / (36.0*pi)
  p0   = 5.0 / (12.0*pi)
  b0   = 1.0 / sqrt(4.0*pi)

  !-- Conditions initiales : (rho, vx, vy, vz, Bx, By, Bz, p, psi)
  do j = 1, ny
    do i = 1, nx
      v(1,i,j) = rho0
      v(2,i,j) = -sin(2.0*pi*y(j))
      v(3,i,j) =  sin(2.0*pi*x(i))
      v(4,i,j) = 0.0
      v(5,i,j) = -b0 * sin(2.0*pi*y(j))
      v(6,i,j) =  b0 * sin(4.0*pi*x(i))
      v(7,i,j) = 0.0
      v(8,i,j) = p0
      v(9,i,j) = 0.0
      call prim_to_cons(v(:,i,j), u(:,i,j))
    end do
  end do

  t = 0.0
  nstep = 0
  do while (t < tmax)
    u0 = u
    call evolve(0.0, 1.0, .true.)
    call evolve(0.5, 0.5, .false.)
    t = t + dt
    nstep = nstep + 1

    !-- amortissement de psi
    do j = 1, ny
      do i = 1, nx
        u(9,i,j) = u(9,i,j) * exp(-alpha_glm * ch * dt / min(dx,dy))
        v(9,i,j) = u(9,i,j)
      end do
    end do

    if (mod(nstep, 50) == 0) print *, 'pas', nstep, ' t =', t, ' max|divB| =', max_divb()
  end do
  print *, 'Fin : t =', t, ' max|divB| =', max_divb()

  open(10, file = 'orszag.txt')
  do j = 1, ny
    do i = 1, nx
      write(10, '(9(es15.5e3,1x))') x(i), y(j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                  v(5,i,j), v(6,i,j), v(8,i,j), v(9,i,j)
    end do
  end do
  close(10)

contains

!-------------------------------------

function max_divb() result(dmax)
  real :: dmax, d
  integer :: i, j
  dmax = 0.0
  do j = 1, ny
    do i = 1, nx
      d = (v(5,i+1,j) - v(5,i-1,j))/(2.0*dx) + (v(6,i,j+1) - v(6,i,j-1))/(2.0*dy)
      dmax = max(dmax, abs(d))
    end do
  end do
end function max_divb

!-------------------------------------

subroutine evolve(c1, c2, set_dt)
  real, intent(in) :: c1, c2
  logical, intent(in) :: set_dt
  real :: vl(nvar), vr(nvar), a, amaxx, amaxy, vm
  real :: fx(nvar, 0:nx, ny), fy(nvar, nx, 0:ny)
  integer :: i, j, k

  !-- Conditions aux limites periodiques
  do j = 1, ny
    v(:, 1-ngc:0,     j) = v(:, nx-ngc+1:nx, j)
    v(:, nx+1:nx+ngc, j) = v(:, 1:ngc,       j)
  end do
  do i = 1-ngc, nx+ngc
    v(:, i, 1-ngc:0    ) = v(:, i, ny-ngc+1:ny)
    v(:, i, ny+1:ny+ngc) = v(:, i, 1:ngc      )
  end do

  !-- Vitesse GLM : la plus grande vitesse d'onde rapide du domaine
  if (set_dt) then
    ch = 0.0
    do j = 1, ny
      do i = 1, nx
        vm = max( abs(v(2,i,j)) + fastspeed(v(:,i,j), 1), &
                  abs(v(3,i,j)) + fastspeed(v(:,i,j), 2) )
        ch = max(ch, vm)
      end do
    end do
  end if

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

end program orszag
