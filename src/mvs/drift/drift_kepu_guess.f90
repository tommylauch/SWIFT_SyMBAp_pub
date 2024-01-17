!*************************************************************************
!                        DRIFT_KEPU_GUESS.F
!*************************************************************************
! Initial guess for solving kepler's equation using universal variables.
!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  initial guess for the value of 
!                                    universal variable
! Author:  Hal Levison & Martin Duncan 
! Date:    3/12/93
! Last revision: April 6/93

subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)
use swift_mod
use orbel_interface
use mvs_interface, except_this_one => drift_kepu_guess
implicit none

real(rk), intent(in)    :: dt,r0,mu,alpha,u

real(rk), intent(inout) :: s

integer(ik)             :: iflg
real(rk)                :: y,sy,cy,sigma,es
real(rk)                :: x,a
real(rk)                :: en,ec,e

!...  Executable code 

   if (alpha.gt.0.0) then
   !...       find initial guess for elliptic motion
      if (dt/r0 .le. 0.4) then
         s = dt/r0-(dt**2*u)/(2.0*r0**3)
         return
      else
         a = mu/alpha
         en = sqrt(mu/(a**3))
         ec = 1.0_rk - r0/a
         es = u/(en*a**2)
         e = sqrt(ec**2 + es**2)
         y = en*dt-es
         call orbel_scget(y,sy,cy)
         sigma = sign(1.0_rk,(es*cy+ec*sy))
         x = y+sigma*0.85_rk*e
         s = x/sqrt(alpha)
      endif
   else
!...       find initial guess for hyperbolic motion.
      call drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
      if (iflg.ne.0) then
         s = dt/r0
      endif
   endif

return
end subroutine drift_kepu_guess
