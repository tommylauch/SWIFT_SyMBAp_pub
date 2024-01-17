!*************************************************************************
!                        DRIFT_KEPU_LAG.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using LAGUERRE'S METHOD
!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflgn          ==>  =0 if converged; !=0 if not
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93

subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
use swift_mod
use mvs_interface, except_this_one => drift_kepu_lag
implicit none

real(rk), intent(in)     :: dt,r0,mu,alpha,u

real(rk), intent(inout)  :: s

real(rk), intent(out)    :: fp,c1,c2,c3
integer(ik), intent(out) :: iflg

integer(ik)              :: nc,ncmax
integer(ik), parameter   :: NTMP = NLAG2+1
real(rk)                 :: ln
real(rk)                 :: x,fpp,ds,c0,f
real(rk)                 :: fdt

!...  Executable code 

!...    To get close approch needed to take lots of iterations if alpha<0
   if(alpha.lt.0.0) then
      ncmax = NLAG2
   else
      ncmax = NLAG2
   endif

   ln = 5.0_rk
!...    start laguere's method
   do nc=0,ncmax
      x = s**2*alpha
      call drift_kepu_stumpff(x,c0,c1,c2,c3)
      c1 = c1*s
      c2 = c2*s**2
      c3 = c3*s**3
      f = r0*c1 + u*c2 + mu*c3 - dt
      fp = r0*c0 + u*c1 + mu*c2
      fpp = (-40.0*alpha + mu)*c1 + u*c0
      ds = - ln*f/(fp + sign(1.d0,fp)*sqrt(abs((ln - 1.0)*             &
                   (ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
      s = s + ds
      fdt = f/dt
!..        quartic convergence
      if (fdt**2.lt.DANBYB**2) then 
         iflg = 0_ik
         return
      endif
!...      Laguerre's method succeeded
   enddo

   iflg = 2_ik

return
end subroutine drift_kepu_lag
