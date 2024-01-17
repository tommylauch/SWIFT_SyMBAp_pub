!*************************************************************************
!                         DRIFT_KEPU_NEW.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using NEWTON'S METHOD
!              Input:
!                  s             ==>  inital value of universal variable
!                  dt            ==>  time step (real scalor)
!                  r0            ==>  Distance between `Sun' and paritcle
!                                      (real scalor)
!                  mu            ==>  Reduced mass of system (real scalor)
!                  alpha         ==>  energy (real scalor)
!                  u             ==>  angular momentun  (real scalor)
!              Output:
!                  s             ==>  final value of universal variable
!                  fp            ==>  f' from p170  
!                                        (real scalors)
!                  c1,c2,c3      ==>  c's from p171-172
!                                        (real scalors)
!                  iflgn          ==>  =0 if converged; !=0 if not
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93

subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)
use swift_mod
use mvs_interface, except_this_one => drift_kepu_new
implicit none

real(rk), intent(in)     :: dt,r0,mu,alpha,u

real(rk), intent(inout)  :: s

real(rk), intent(out)    :: fp,c1,c2,c3
integer(ik), intent(out) :: iflgn

integer(ik)              :: nc
real(rk)                 :: x,c0,ds
real(rk)                 :: f,fpp,fppp,fdt

!....  Executable code 

   do nc=0,6
      x = s**2*alpha
      call drift_kepu_stumpff(x,c0,c1,c2,c3)
      c1 = c1*s 
      c2 = c2*s**2
      c3 = c3*s**3
      f = r0*c1 + u*c2 + mu*c3 - dt
      fp = r0*c0 + u*c1 + mu*c2
      fpp = (-r0*alpha + mu)*c1 + u*c0
      fppp = (- r0*alpha + mu)*c0 - u*alpha*c1
      ds = - f/fp
      ds = - f/(fp + ds*fpp*0.5_rk)
      ds = -f/(fp + ds*fpp*0.5_rk + ds*ds*fppp*0.5_rk*ONETHRD)
      s = s + ds
      fdt = f/dt

      !...      quartic convergence
      if (fdt**2.lt.DANBYB**2) then 
         iflgn = 0_ik
         return
      endif
   !...     newton's method succeeded

   enddo

!...     newton's method failed
   iflgn = 1_ik

return
end subroutine drift_kepu_new
