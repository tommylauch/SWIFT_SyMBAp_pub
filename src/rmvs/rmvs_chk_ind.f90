!*************************************************************************
!                            RMVS_CHK_IND.F
!*************************************************************************
!  Subroutine to check if a test particle and planet
!  are having or **will** have an encounter 
!  in the next timestep. 
!             Input:
!                 xr           ==>  relative position of tp wrt planet
!                                   (real scalar)
!                 vxr          ==>  relative velocity of tp wrt planet
!                                   (real scalor)
!                 dt           ==>  time step (real scalor)
!                 r2crit       ==> boundary of outer enc region
!                                   (real scalor)
!                 r2critp      ==> boundary of inner (planocentric) enc region
!                                   (real scalor)
!             Output:
!                 iflag        ==> encounter?  =  0 no
!                                              =  1 yes, in outer region
!                                              = -1 yes, in inner region
! Remarks: Based on Hal's wiscl_fk.f' but origonaly written by Martin Duncan
! Authors:  Hal Levison 
! Date:    2/19/93
! Last revision: 

subroutine rmvs_chk_ind(xr,vxr,dt,r2crit,r2critp,iflag)
use swift_mod
implicit none

real(rk), intent(in)     :: xr(:),vxr(:),dt,r2crit,r2critp

integer(ik), intent(out) :: iflag

real(rk)                 :: r2,v2,vdotr,tmin,r2min

!...  Executable code

!...    First check if we're already in the encounter region. If so return
!.             with flag set to one.
   r2 = dot_product(xr(1:3),xr(1:3))
   if (r2 .le. r2critp) then
      iflag = -1_ik
      return
   endif

!...    If we're heading outward, use r2 to calc iflag
vdotr = dot_product(xr(1:3),vxr(1:3))
   if (vdotr.gt.0.0_rk) then
      if (r2.ge.r2crit) then
         iflag = 0_ik
      else
         iflag = 1_ik
      endif
      return
   endif      

!...    We're not yet inside and are converging so we need to calc. the
!.           minimum separation attained in time dt.
   v2 = dot_product(vxr(1:3),vxr(1:3))
   tmin = -vdotr/v2

   if (tmin.lt.dt) then
      r2min = r2-(vdotr**2)/v2
   else
      r2min = r2+2.0_rk*vdotr*dt + v2*dt**2
   endif

   r2min = min(r2min,r2)     ! really make sure

   if (r2min.le.r2critp) then
      iflag = -1_ik
   else if (r2min.le.r2crit) then
      iflag = 1_ik
   else 
      iflag = 0_ik
   endif

return
end subroutine rmvs_chk_ind
