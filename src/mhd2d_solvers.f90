module mhd2d_solvers

  implicit none

  real, parameter :: gamma = 5.0/3.0
  integer, parameter :: nvar = 9
  !-- 1=rho, 2=vx, 3=vy, 4=vz, 5=Bx, 6=By, 7=Bz, 8=p, 9=psi (GLM)

  real :: ch = 1.0     !-- vitesse de nettoyage GLM (fixee a chaque pas par le programme)

contains

!-------------------------------------
function reconstruction(a1, a2, a3) result(rec)
  real, intent(in) :: a1, a2, a3
  real :: rec, r, phi
  real, parameter :: eps = 1e-10

  r = (a2 - a1) / (a3 - a2 + eps)
  phi = max(0.0, min(1.0, r))          ! minmod
  rec = a2 + 0.5 * phi * (a3 - a2)
end function reconstruction

!-------------------------------------
function fastspeed(v, dir) result(cf)
  !-- vitesse magnetosonique rapide selon la direction dir (1=x, 2=y)
  real, intent(in) :: v(nvar)
  integer, intent(in) :: dir
  real :: cf, a2, b2, bn2, disc

  a2 = gamma * v(8) / v(1)
  b2 = (v(5)**2 + v(6)**2 + v(7)**2) / v(1)
  if (dir == 1) then
    bn2 = v(5)**2 / v(1)
  else
    bn2 = v(6)**2 / v(1)
  end if
  disc = max( (a2 + b2)**2 - 4.0*a2*bn2, 0.0 )
  cf = sqrt( 0.5 * ( (a2 + b2) + sqrt(disc) ) )
end function fastspeed

!-------------------------------------
subroutine llf(vl, vr, dir, f, a)
  real, intent( in) :: vl(nvar), vr(nvar)
  integer, intent(in) :: dir
  real, intent(out) :: f(nvar), a
  real :: fl(nvar), fr(nvar), ul(nvar), ur(nvar)
  real :: lpl, lml, lpr, lmr

  call prim_to_cons(vl, ul)
  call prim_to_cons(vr, ur)
  call flux(vl, dir, fl, lpl, lml)
  call flux(vr, dir, fr, lpr, lmr)

  a = max( abs(lpl), abs(lml), abs(lpr), abs(lmr), ch )   ! ch : ondes GLM
  f = 0.5 * (fl + fr) - 0.5 * a * (ur - ul)
end subroutine llf

!-------------------------------------
subroutine prim_to_cons(v, u)
  real, intent( in) :: v(nvar)
  real, intent(out) :: u(nvar)
  real :: b2

  b2 = v(5)**2 + v(6)**2 + v(7)**2
  u(1) = v(1)
  u(2) = v(1) * v(2)
  u(3) = v(1) * v(3)
  u(4) = v(1) * v(4)
  u(5) = v(5)
  u(6) = v(6)
  u(7) = v(7)
  u(8) = v(8)/(gamma-1.0) + 0.5*v(1)*(v(2)**2 + v(3)**2 + v(4)**2) + 0.5*b2
  u(9) = v(9)
end subroutine prim_to_cons

!-------------------------------------
subroutine cons_to_prim(u, v)
  real, intent( in) :: u(nvar)
  real, intent(out) :: v(nvar)
  real :: b2

  b2 = u(5)**2 + u(6)**2 + u(7)**2
  v(1) = u(1)
  v(2) = u(2)/u(1)
  v(3) = u(3)/u(1)
  v(4) = u(4)/u(1)
  v(5) = u(5)
  v(6) = u(6)
  v(7) = u(7)
  v(8) = (gamma-1.0) * ( u(8) - 0.5*u(1)*(v(2)**2 + v(3)**2 + v(4)**2) - 0.5*b2 )
  v(9) = u(9)
end subroutine cons_to_prim

!-------------------------------------
subroutine flux(v, dir, f, lp, lm)
  !-- flux physique selon x (dir=1) ou y (dir=2), avec les termes GLM
  real, intent( in) :: v(nvar)
  integer, intent(in) :: dir
  real, intent(out) :: f(nvar), lp, lm
  real :: rho, vx, vy, vz, bx, by, bz, p, psi
  real :: b2, ptot, e, vdotb, vn, cf

  rho = v(1); vx = v(2); vy = v(3); vz = v(4)
  bx = v(5);  by = v(6); bz = v(7); p = v(8); psi = v(9)
  b2 = bx**2 + by**2 + bz**2
  ptot = p + 0.5*b2
  e = p/(gamma-1.0) + 0.5*rho*(vx**2 + vy**2 + vz**2) + 0.5*b2
  vdotb = vx*bx + vy*by + vz*bz

  if (dir == 1) then
    vn = vx
    f(1) = rho*vx
    f(2) = rho*vx**2 + ptot - bx**2
    f(3) = rho*vx*vy - bx*by
    f(4) = rho*vx*vz - bx*bz
    f(5) = psi                       ! GLM : flux de Bx
    f(6) = by*vx - bx*vy
    f(7) = bz*vx - bx*vz
    f(8) = (e + ptot)*vx - bx*vdotb
    f(9) = ch**2 * bx                ! GLM : flux de psi
  else
    vn = vy
    f(1) = rho*vy
    f(2) = rho*vx*vy - bx*by
    f(3) = rho*vy**2 + ptot - by**2
    f(4) = rho*vy*vz - by*bz
    f(5) = bx*vy - by*vx
    f(6) = psi                       ! GLM : flux de By
    f(7) = bz*vy - by*vz
    f(8) = (e + ptot)*vy - by*vdotb
    f(9) = ch**2 * by
  end if

  cf = fastspeed(v, dir)
  lp = vn + cf
  lm = vn - cf
end subroutine flux

end module mhd2d_solvers
