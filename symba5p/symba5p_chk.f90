!*************************************************************************
!                            SYMBA5P_CHK.F
!*************************************************************************
! This subroutine checks to see if there are encounters
!             Input:
!                 rhill         ==>  Radius of hill sphere (real array)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip1,ip2       ==>  The two bodies to check (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh,yh,zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 dt            ==>  time step  (real scalor)
!                 irec          ==>  current recursion level (int scalar)
!             Output:
!                 icflg         ==> ecounter?  = 1 Yes
!                                              =  0 No (integer scalar)  
!                 svdotr        ==> = .true. if i,j are receding
!                                   = .false is approaching
!                                     (logical*1 scalar)
! Remarks: Based on plh_chk.f.  Same as symba_chk.f
! Authors:  Hal Levison
! Date:   3/20/97
! Last revision: 

subroutine symba5p_chk(rhill,nbod,ip1,ip2,mass,xh,vxh,dt,              &
                       irec,icflg,svdotr)
implicit none
use swift_mod
use sybam5p_mod
use rmvs_interface

integer(ik), intent(in)  :: nbod,irec,ip1,ip2
real(rk), intent(in)     :: mass(:),xh(:,:),dt
real(rk), intent(in)     :: vxh(:,:),rhill(:)

integer(ik), intent(out) :: icflg
logical(ik), intent(out) :: svdotr

real(rk)                 :: r2crit,r2critp,rcrit
real(rk)                 :: xr(3),vxr(3)
real(rk)                 :: vdotr

!...  Executable code 

   rcrit = (rhill(ip1)+rhill(ip2)) * RHSCALE * (RSHELL**(irec))
   r2crit = rcrit*rcrit
   r2critp = -1.0_rk          ! not used here

   xr(:) = xh(:,ip2) - xh(:,ip1)
   vxr(:) = vxh(:,ip2) - vxh(:,ip1)
   call rmvs_chk_ind(xr,vxr,dt,r2crit,r2critp,icflg)

   vdotr = dot_product(xr,vxr)
   svdotr = (vdotr.lt.0.0_rk)

return
end subroutine symba5p_chk


