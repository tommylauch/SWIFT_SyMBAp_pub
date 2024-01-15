c*************************************************************************
c                            UTIL_HILLS.F
c*************************************************************************
c This subroutine calculates the hill's sphere for the planets
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  initial position in helio coord 
c                                    (real array)
c                 vxh           ==>  initial velocity in helio coord 
c                                    (real array)
c             Output:
c                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
c                                    (real array)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 1/6/97

subroutine util_hills(nbod,mass,xh,vxh,r2hill)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)

real(rk), intent(out)   :: r2hill(:)

integer(ik)             :: i
real(rk)                :: mu,energy,ap,rhil,r,v2

c...  Executable code 

do i=2,nbod
   if (mass(i) .ne. 0.0_rk) then
      mu = mass(1)*mass(i)/(mass(1)+mass(i))
      r = sqrt(dot_product(xh,xh))
      v2 = dot_product(vxh,vxh)
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
