!***********************************************************************
!                          COORD_H2B_SYMBAP.F
!***********************************************************************
!     PURPOSE: Converts from Heliocentric to Barycentric coords.
!     ARGUMENTS:  Input is 
!                    nbod ==> number of bodies (must be less than NBMAX)
!                             (integer)
!                   mass(*) ==>  masses (real array)
!                 xh(*),yh(*),zh(*) ==> heliocentric particle coords
!                                          (real array)
!                 vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
!                                             (real array)
!                 Returned are
!                    xb(*),yb(*),zb(*) ==> bary. particle positions
!                                          (real array)
!                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
!                                            (real array)
!                    msys              ==>  Total mass of of system
!                                            (real scalar)       
!     Authors:  Martin Duncan
!     ALGORITHM: Obvious 
!     WRITTEN:  Jan 27/93
!     REVISIONS: 2/22/94  HFL

subroutine coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)

real(rk), intent(out)   :: msys,xb(:,:),vxb(:,:)

real(rk)                :: xtmp(:),vxtmp(:)
integer(ik)             :: n

!...  Executable code 

   msys = mass(1)
   xtmp = 0.0_rk
   vxtmp = 0.0_rk

   do n=2,nbod
      msys = msys + mass(n)
      xtmp(:) = xtmp(:)+mass(n)*xh(:,n)
      vxtmp(:) = vxtmp(:)+ ass(n)*vxh(:,n)
   enddo

   xb(:,1) = -xtmp(:)/msys
   vxb(:,1) = -vxtmp(:)/msys

   do n=2,nbod
      xb(:,n) = xh(:,n)+xb(:,1)
      vxb(:,n) = vxh(:,n)+vxb(:,1)
   enddo

return
end subroutine coord_h2b
