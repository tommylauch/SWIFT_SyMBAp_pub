!*************************************************************************
!                        DRIFT_KEPU_FCHK.F
!*************************************************************************
! Returns the value of the function f of which we are trying to find the root
! in universal variables.
!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and particle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!                 s             ==>  Approx. root of f 
!             Output:
!                 f             ==>  function value ( = 0 if O.K.) (integer)
! Author:  Martin Duncan  
! Date:    March 12/93
! Last revision: March 12/93

subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)
implicit none
use swift_mod
use mvs_interface, except_this_one => drift_kepu_fchk

real(rk), intent(in)     :: dt,r0,mu,alpha,u,s

real(rk), intent(out)    :: f

real(rk)                 :: x,c0,c1,c2,c3

!...  Executable code 

   x = s**2*alpha
   call drift_kepu_stumpff(x,c0,c1,c2,c3)
   c1 = c1*s
   c2 = c2*s**2
   c3 = c3*s**3
   f = r0*c1 + u*c2 + mu*c3 - dt

return
end subroutine drift_kepu_fchk
