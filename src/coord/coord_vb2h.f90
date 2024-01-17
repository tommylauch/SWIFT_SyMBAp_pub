!***********************************************************************
!                          COORD_VB2H.F
!***********************************************************************
!     PURPOSE: Converts from Barycentric to Helio coords.
!               Velocity only
!     ARGUMENTS:  Input is 
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
!                                 SYMMETRY IN SUBROUTINE CALLS
!                 vxb(*) ==> Barycentric particle velocities
!                                             (real array)
!                 Returned are
!                    vxh(*) ==> Helio particle velocities
!                                            (real array)
!       
!     ALGORITHM: Obvious 
!     Authors:  Hal Levison
!     WRITTEN:  11/14/96
!     REVISIONS: 11/21/96

subroutine coord_vb2h(nbod,mass,vxb,vxh)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:)

real(rk), intent(inout) :: vxb(:,:)

real(rk), intent(out)   :: vxh(:,:)

integer(ik)             :: i

!...  Executable code

   vxb(:,1) = -mass(2)*vxb(:,2)

   do i=3,nbod
      vxb(:,1) = vxb(:,1)-mass(i)*vxb(:,i)
   enddo

   vxb(:,1) = vxb(:,1)/mass(1)

   do i=2,nbod
      vxh(:,i) = vxb(:,i)-vxb(:,1)
   enddo

return
end subroutine coord_vb2h
