!***********************************************************************
!                          COORD_VH2B.F
!***********************************************************************
!     PURPOSE: Converts from Heliocentric to Barycentric coords. 
!              Velocity only
!     ARGUMENTS:  Input is 
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                 vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
!                                             (real array)
!                 Returned are
!                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
!                                            (real array)
!                    msys              ==>  Total mass of of system
!                                            (real scalar)       
!     Authors:  Hal Levison
!     ALGORITHM: Obvious 
!     WRITTEN:  11/14/96
!     REVISIONS: 11/21/96

subroutine coord_vh2b(nbod,mass,vxh,vxb)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),vxh(:,:)

real(rk), intent(out)   :: vxb(:,:)

real(rk)                :: msys,vxtmp(3)
integer(ik)             :: n

!...  Executable code

   msys = mass(1)
   vxtmp = 0.0_rk

   do n=2,nbod
      msys = msys+mass(n)
      vxtmp(:) = vxtmp(:) + mass(n)*vxh(:,n)
   enddo

   vxb(:,1) = -vxtmp(:)/msys

   do n=2,nbod
     vxb(:,n) = vxh(:,n) + vxb(:,1)
   enddo

return
end subroutine coord_vh2b
