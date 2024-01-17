!*************************************************************************
!                            UTIL_HILLS1.F
!*************************************************************************
! This subroutine calculates the hill's sphere for the planets
!             Input:
!                 msun          ==>  mass of sun (real scalar)
!                 mpl           ==>  mass of sun (real scalar)
!                 xh            ==>  position of pl in helio coord 
!                                    (real scalars)
!                 vxh           ==>  velocity of pl in helio coord 
!                                    (real scalars)
!             Output:
!                  rhill        ==>  the radius of planet's hill's sphere 
!                                    (real scalar)
! Remarks: Based on util_hill
! Authors:  Hal Levison 
! Date:    1/8/97
! Last revision: 

subroutine util_hills1(msun,mpl,xh,vxh,rhill) 
use swift_mod
implicit none

real(rk), intent(in)       :: msun,mpl,xh(:),vxh(:)

real(rk), intent(out)      :: rhill

real(rk)                   :: mu,energy,ap,r,v2

!...  Executable code 

   mu = msun*mpl/(msun+mpl)
   r = sqrt(dot_product(xh(1:3),xh(1:3)))
   v2 = dot_product(vxh(1:3),vxh(1:3))
   energy = -1.0_rk*msun*mpl/r+0.5_rk*mu*v2
   ap = -1.0_rk*msun*mpl/(2.0_rk*energy)
   rhill = ap*(((mu/msun)/3.0_rk)**ONETHRD)

return
end subroutine util_hills1
