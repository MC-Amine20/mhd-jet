program jetmhd
  !-- Jet MHD 2D (tranche) avec champ axial uniforme, nettoyage GLM de div B
  !-- Compilation : gfortran -O2 -fdefault-real-8 mhd2d_solvers.f90 jetmhd.f90 -o jetmhd
  !-- Execution   : .\jetmhd 1.0     (beta du plasma = 1)
  !--               .\jetmhd 0       (aucun champ magnetique : run de controle)

  use mhd2d_solvers
  implicit none

  integer, parameter :: nx = 240, ny = 160, ngc = 2
  real, parameter :: cfl = 0.4
  real, parameter :: xmax = 12.0, ymax = 8.0
  real, parameter :: tmax = 2.4
  integer, parameter :: nout = 24
  real, parameter :: alpha_glm = 0.1

  !-- Milieu ambiant et jet
  real, parameter :: rho_a = 1.0, p_a = 1.0
  real, parameter :: eta = 0.1, mach = 5.0, rj = 1.0
  real, parameter :: pfloor = 1.0e-6, rhofloor = 1.0e-6

  real :: x(nx), y(ny)
  real :: v(nvar, 1-ngc:nx+ngc, 1-ngc:ny+ngc)
  real :: u(nvar, nx, ny), u0(nvar, nx, ny)
  real :: vjet(nvar), vamb(nvar)
  real :: t, dt, dx, dy, yc, rho_j, p_j, vj, dtout, tnext
  real :: beta, b0
  character(len=32) :: arg, tag
  integer :: i, j, n, nstep

  !-- beta du plasma lu sur la ligne de commande (defaut 1.0 ; 0 = sans champ)
  beta = 1.0
  tag = '1.0'
  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) then
    read(arg, *) beta
    tag = trim(arg)
  end if
  if (beta <= 0.0) then
    b0 = 0.0
  else
    b0 = sqrt(2.0 * p_a / beta)          !-- beta = 2 p / B^2
  end if

  dx = xmax / real(nx)
  dy = ymax / real(ny)
  yc = 0.5 * ymax
  do i = 1, nx
    x(i) = (real(i) - 0.5) * dx
  end do
  do j = 1, ny
    y(j) = (real(j) - 0.5) * dy
  end do

  !-- Etats : (rho, vx, vy, vz, Bx, By, Bz, p, psi)
  !-- Champ axial uniforme Bx = b0 dans le jet et l'ambiant (meme pression totale)
  rho_j = eta * rho_a
  p_j   = p_a
  vj    = mach * sqrt(gamma * p_j / rho_j)
  vjet  = [ rho_j, vj,  0.0, 0.0, b0, 0.0, 0.0, p_j, 0.0 ]
  vamb  = [ rho_a, 0.0, 0.0, 0.0, b0, 0.0, 0.0, p_a, 0.0 ]
  print *, 'beta =', beta, '  B0 =', b0, '  vitesse du jet =', vj
  if (b0 > 0.0) print *, 'Mach alfvenique du jet =', vj / (b0/sqrt(rho_j))

  do j = 1, ny
    do i = 1, nx
      v(:,i,j) = vamb
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

    !-- amortissement de psi
    do j = 1, ny
      do i = 1, nx
        u(9,i,j) = u(9,i,j) * exp(-alpha_glm * ch * dt / min(dx,dy))
        v(9,i,j) = u(9,i,j)
      end do
    end do

    if (mod(nstep, 200) == 0) print *, 'pas', nstep, ' t =', t
    if (t >= tnext - 1.0e-10) then
      n = n + 1
      call write_snapshot(n)
      tnext = tnext + dtout
    end if
  end do
  print *, 'Termine : ', nstep, ' pas'

contains

!-------------------------------------

subroutine write_snapshot(num)
  integer, intent(in) :: num
  character(len=64) :: fname
  integer :: i, j

  write(fname, '(a,a,a,i3.3,a)') 'jetmhd_', trim(tag), '_', num, '.txt'
  open(10, file = trim(fname))
  do i = 1, nx
    do j = 1, ny
      write(10, '(8(es15.5e3,1x))') x(i), y(j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                    v(8,i,j), v(5,i,j), v(6,i,j)
    end do
  end do
  close(10)
end subroutine write_snapshot

!-------------------------------------

subroutine evolve(c1, c2, set_dt)
  real, intent(in) :: c1, c2
  logical, intent(in) :: set_dt
  real :: vl(nvar), vr(nvar), a, amaxx, amaxy, vm
  real :: fx(nvar, 0:nx, ny), fy(nvar, nx, 0:ny)
  integer :: i, j, k

  !-- Bord gauche (x=0) : injection du jet ou mur (vx inverse, Bx impose)
  do j = 1, ny
    if (abs(y(j) - yc) <= rj) then
      v(:, 0,  j) = vjet
      v(:, -1, j) = vjet
    else
      v(:, 0,  j) = v(:, 1, j);  v(2, 0,  j) = -v(2, 1, j);  v(5, 0,  j) = b0
      v(:, -1, j) = v(:, 2, j);  v(2, -1, j) = -v(2, 2, j);  v(5, -1, j) = b0
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
      if (v(1,i,j) < rhofloor .or. v(8,i,j) < pfloor) then
        v(1,i,j) = max(v(1,i,j), rhofloor)
        v(8,i,j) = max(v(8,i,j), pfloor)
        call prim_to_cons(v(:,i,j), u(:,i,j))
      end if
    end do
  end do

end subroutine evolve

end program jetmhd
