module mhd_solvers

  implicit none

  real, parameter :: gamma = 2.0
  real, parameter :: bx = 0.75      !-- Bx constant en 1D
  integer, parameter :: nvar = 7

contains

!-------------------------------------
function reconstruction(a1, a2, a3) result(rec)
  real, intent(in) :: a1, a2, a3
  real :: rec, r, phi
  real, parameter :: eps = 1e-10

  r = (a2 - a1) / (a3 - a2 + eps)
  phi = max(0.0, min(1.0, r))          ! minmod (plus robuste en MHD)
  rec = a2 + 0.5 * phi * (a3 - a2)
end function reconstruction

!-------------------------------------
subroutine llf(vl, vr, f, a)
  real, intent( in) :: vl(nvar), vr(nvar)
  real, intent(out) :: f(nvar), a
  real :: fl(nvar), fr(nvar), ul(nvar), ur(nvar)
  real :: lpl, lml, lpr, lmr

  call prim_to_cons(vl, ul)
  call prim_to_cons(vr, ur)
  call flux(vl, fl, lpl, lml)
  call flux(vr, fr, lpr, lmr)

  a = max( abs(lpl), abs(lml), abs(lpr), abs(lmr) )
  f = 0.5 * (fl + fr) - 0.5 * a * (ur - ul)
end subroutine llf

!-------------------------------------
subroutine prim_to_cons(v, u)
  real, intent( in) :: v(nvar)
  real, intent(out) :: u(nvar)
  real :: rho, vx, vy, vz, by, bz, p, b2

  rho = v(1); vx = v(2); vy = v(3); vz = v(4)
  by = v(5);  bz = v(6); p = v(7)
  b2 = bx**2 + by**2 + bz**2

  u(1) = rho
  u(2) = rho * vx
  u(3) = rho * vy
  u(4) = rho * vz
  u(5) = by
  u(6) = bz
  u(7) = p/(gamma-1.0) + 0.5*rho*(vx**2 + vy**2 + vz**2) + 0.5*b2
end subroutine prim_to_cons

!-------------------------------------
subroutine cons_to_prim(u, v)
  real, intent( in) :: u(nvar)
  real, intent(out) :: v(nvar)
  real :: rho, vx, vy, vz, by, bz, b2

  rho = u(1)
  vx = u(2)/rho; vy = u(3)/rho; vz = u(4)/rho
  by = u(5);     bz = u(6)
  b2 = bx**2 + by**2 + bz**2

  v(1) = rho
  v(2) = vx
  v(3) = vy
  v(4) = vz
  v(5) = by
  v(6) = bz
  v(7) = (gamma-1.0) * ( u(7) - 0.5*rho*(vx**2+vy**2+vz**2) - 0.5*b2 )
end subroutine cons_to_prim

!-------------------------------------
subroutine flux(v, f, lp, lm)
  real, intent( in) :: v(nvar)
  real, intent(out) :: f(nvar), lp, lm
  real :: rho, vx, vy, vz, by, bz, p, b2, e, vdotb
  real :: a2, bb2, bbx2, cf

  rho = v(1); vx = v(2); vy = v(3); vz = v(4)
  by = v(5);  bz = v(6); p = v(7)
  b2 = bx**2 + by**2 + bz**2
  e  = p/(gamma-1.0) + 0.5*rho*(vx**2+vy**2+vz**2) + 0.5*b2
  vdotb = vx*bx + vy*by + vz*bz

  f(1) = rho*vx
  f(2) = rho*vx**2 + p + 0.5*b2 - bx**2
  f(3) = rho*vx*vy - bx*by
  f(4) = rho*vx*vz - bx*bz
  f(5) = by*vx - bx*vy
  f(6) = bz*vx - bx*vz
  f(7) = (e + p + 0.5*b2)*vx - bx*vdotb

  !-- vitesse magnetosonique rapide
  a2   = gamma*p/rho
  bb2  = b2/rho
  bbx2 = bx**2/rho
  cf = sqrt( 0.5*( (a2+bb2) + sqrt( (a2+bb2)**2 - 4.0*a2*bbx2 ) ) )

  lp = vx + cf
  lm = vx - cf
end subroutine flux

end module mhd_solvers