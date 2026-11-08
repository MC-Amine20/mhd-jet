module hydro2d_solvers

  implicit none

  real, parameter :: gamma = 1.4
  integer, parameter :: nvar = 4

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

  a = max( abs(lpl), abs(lml), abs(lpr), abs(lmr) )
  f = 0.5 * (fl + fr) - 0.5 * a * (ur - ul)
end subroutine llf

!-------------------------------------
subroutine prim_to_cons(v, u)
  real, intent( in) :: v(nvar)
  real, intent(out) :: u(nvar)

  u(1) = v(1)
  u(2) = v(1) * v(2)
  u(3) = v(1) * v(3)
  u(4) = v(4)/(gamma-1.0) + 0.5*v(1)*(v(2)**2 + v(3)**2)
end subroutine prim_to_cons

!-------------------------------------
subroutine cons_to_prim(u, v)
  real, intent( in) :: u(nvar)
  real, intent(out) :: v(nvar)

  v(1) = u(1)
  v(2) = u(2)/u(1)
  v(3) = u(3)/u(1)
  v(4) = (gamma-1.0) * ( u(4) - 0.5*u(1)*(v(2)**2 + v(3)**2) )
end subroutine cons_to_prim

!-------------------------------------
subroutine flux(v, dir, f, lp, lm)
  !-- dir = 1 : flux selon x ; dir = 2 : flux selon y
  real, intent( in) :: v(nvar)
  integer, intent(in) :: dir
  real, intent(out) :: f(nvar), lp, lm
  real :: rho, vx, vy, p, vn, cs, e

  rho = v(1); vx = v(2); vy = v(3); p = v(4)
  e  = p/(gamma-1.0) + 0.5*rho*(vx**2 + vy**2)
  cs = sqrt(gamma*p/rho)

  if (dir == 1) then
    vn = vx
  else
    vn = vy
  end if

  f(1) = rho*vn
  f(2) = rho*vx*vn
  f(3) = rho*vy*vn
  f(4) = (e + p)*vn
  if (dir == 1) then
    f(2) = f(2) + p
  else
    f(3) = f(3) + p
  end if

  lp = vn + cs
  lm = vn - cs
end subroutine flux

end module hydro2d_solvers