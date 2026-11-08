module hydro_solvers

  implicit none
  
  real, parameter :: gamma = 1.4

  integer, parameter :: nvar = 3
  !-- 1 = rh, 2 = vx, 3 = pr
  
contains

!-------------------------------------

function reconstruction(a1, a2, a3) result(rec)
!-- Reconstruction of v_{i+1/2} given a1=v_{i-1}, a2=v_{i}, a3=v_{i+1}

  real, intent(in) :: a1, a2, a3
  real :: rec, r, phi
  real, parameter :: eps = 1e-10


  r = (a2 - a1) / (a3 - a2 + eps)

  ! Van Leer limiter
  phi = (r + abs(r)) / (1.0 + abs(r))

  ! Linear reconstruction TVD
  rec = a2 + 0.5 * phi * (a3 - a2)

  ! Alternative plus dissipative (minmod) – décommente si besoin
  ! phi = max(0.0, min(1.0, r))
  ! rec = a2 + 0.5 * phi * (a3 - a2)
 
  
end function reconstruction

!-------------------------------------

subroutine llf(vl, vr, f, a)
!-- Numerical flux LLF
  
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
!-- From primitive v to conserved u

  real, intent( in) :: v(nvar)
  real, intent(out) :: u(nvar)
  real :: rho, vx, P 
  rho = v(1)
  vx = v(2)
  P = v(3)
  
  u(1) = rho 
  u(2) = rho * vx 
  u(3) = P/(gamma-1) + 0.5 * rho * vx**2


end subroutine prim_to_cons

!-------------------------------------

subroutine cons_to_prim(u, v)
!-- From conserved u to primitive v

  real, intent( in) :: u(nvar)
  real, intent(out) :: v(nvar)
  real :: rho, rho_vx, E, P, vx 
  
  rho = u(1)
  rho_vx = u(2)
  E = u(3)

  vx = rho_vx/rho
  P = (gamma - 1) * (E - 0.5 * rho * vx**2)

  v(1) = rho
  v(2) = vx
  v(3) = P

  
end subroutine cons_to_prim

!-------------------------------------

subroutine flux(v, f, lp, lm)
!-- From primitive v to (physical) flux f and lambda^+, lambda^-

  real, intent( in) :: v(nvar)
  real, intent(out) :: f(nvar), lp, lm
  
  real :: rho, cs, vx, P, E
  rho = v(1)
  vx = v(2)
  p = v(3)
  
  cs = sqrt(gamma*P/rho) ! sound speed
  E = P/(gamma-1) + 0.5 * rho * vx**2 !energy

  !flux values
  f(1) = rho * vx
  f(2) = rho * vx**2 + P
  f(3) = vx*E + vx*P 

  !eigen value
  lp = cs + vx
  lm = vx - cs

end subroutine flux

!-------------------------------------
  
end module hydro_solvers
