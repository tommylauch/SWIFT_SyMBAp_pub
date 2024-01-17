!*************************************************************************
!                            UTIL_HILLS.F
!*************************************************************************
! This subroutine calculates the hill's sphere for the planets
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>  initial position in helio coord 
!                                    (real array)
!                 vxh           ==>  initial velocity in helio coord 
!                                    (real array)
!             Output:
!                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
!                                    (real array)
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/19/93
! Last revision: 1/6/97

subroutine util_hills(nbod,mass,xh,vxh,r2hill)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)

real(rk), intent(out)   :: r2hill(:)

integer(ik)             :: i
real(rk)                :: mu,energy,ap,rhil,r,v2

!...  Executable code 

   do i=2,nbod
      if (mass(i).ne.0.0_rk) then
         mu = mass(1)*mass(i)/(mass(1)+mass(i))
         r = sqrt(dot_product(xh(1:3,i),xh(1:3,i)))
         v2 = dot_product(vxh(1:3,i),vxh(1:3,i))
         energy = -1.0_rk*mass(1)*mass(i)/r + 0.5_rk*mu*v2
         ap = -1.0_rk*mass(1)*mass(i)/(2.0_rk*energy)
         rhil = ap*(((mu/mass(1))/3.0_rk)**ONETHRD)
         r2hill(i) = rhil**2
      else
         r2hill(i) = 0.0_rk
      endif
   enddo

   r2hill(1) = 0.0_rk

return
end subroutine util_hills
