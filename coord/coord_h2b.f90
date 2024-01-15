c***********************************************************************
c                          COORD_H2B_SYMBAP.F
c***********************************************************************
c     PURPOSE: Converts from Heliocentric to Barycentric coords.
c     ARGUMENTS:  Input is 
c                    nbod ==> number of bodies (must be less than NBMAX)
c                             (integer)
c                   mass(*) ==>  masses (real array)
c                 xh(*),yh(*),zh(*) ==> heliocentric particle coords
c                                          (real array)
c                 vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
c                                             (real array)
c                 Returned are
c                    xb(*),yb(*),zb(*) ==> bary. particle positions
c                                          (real array)
c                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
c                                            (real array)
c                    msys              ==>  Total mass of of system
c                                            (real scalar)       
c     Authors:  Martin Duncan
c     ALGORITHM: Obvious 
c     WRITTEN:  Jan 27/93
c     REVISIONS: 2/22/94  HFL

subroutine coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)

real(rk), intent(out)   :: xb(:,:),vxb(:,:)

real(rk)                :: msys,xtmp(:),vxtmp(:)
integer(ik)             :: n

c...  Executable code 

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
