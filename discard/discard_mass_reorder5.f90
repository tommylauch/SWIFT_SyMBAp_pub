!*************************************************************************
!                            DISCARD_MASS_REORDER5.F
!*************************************************************************
! Remove a massive body
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ip            ==>  planets to remove (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>   position in helio coord 
!                                    (real arrays)
!                 vxh           ==>   pl vel in helio coord 
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isperih       ==> heliocentric peri flags. (real array)
!             Output:
!                 ip            ==>  planets to remove (int scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>   position in helio coord 
!                                    (real arrays)
!                 vxh           ==>   pl vel in helio coord 
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isperih       ==> heliocentric peri flags. (real array)
! Remarks: 
! Authors:  Hal Levison 
! Date:    1/2/97
! Last revision: 5/13/99

subroutine discard_mass_reorder5(ip,nbod,mass,xh,vxh,rpl,rhill,isperih)
implicit none
use swift_mod

integer(ik), intent(in)    :: ip

integer(ik), intent(inout) :: nbod
real(rk), intent(inout)    :: mass(:),xh(:,:)
real(rk), intent(inout)    :: vxh(:,:),rpl(:)
real(rk), intent(inout)    :: rhill(:)
integer(ik), intent(inout) :: isperih(:)

integer(ik)                :: i,j

!...  Executable code 

   do i=ip,nbod-1
      xh(:,i) = xh(:,i+1)
      vxh(:,i) = vxh(:,i+1)
      mass(i) = mass(i+1)
      rpl(i) = rpl(i+1)
      rhill(i) = rhill(i+1)
      isperih(i) = isperih(i+1)
   enddo
   nbod = nbod-1

return
end subroutine discard_mass_reorder5
