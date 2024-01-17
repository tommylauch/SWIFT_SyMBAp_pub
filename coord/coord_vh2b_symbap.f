c***********************************************************************
c                          COORD_VH2B_SYMBAP.F
c***********************************************************************
*     PURPOSE: Converts from Heliocentric to Barycentric coords. 
*              Velocity only
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                   mass(*) ==>  masses (real array)
*                 vxh(*),vyh(*),vzh(*) ==> heliocentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxb(*),vyb(*),vzb(*) ==> bary. particle velocities
*                                            (real array)
*                    msys              ==>  Total mass of of system
*                                            (real scalar)       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  11/14/96
*     REVISIONS: 11/21/96

      subroutine coord_vh2b_symbap(nbod,mass,vxh,vxb,msys)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(*),vxh(3,*)

c...  Outputs:
      real*8 vxb(3,*)

c...  Internals:
      real*8 msys,vxtmp(3)
      integer n

c----
c...  Executable code 

      msys = mass(1)
      vxtmp = 0.d0

      do n=2,nbod
         msys = msys + mass(n)
         vxtmp(:) = vxtmp(:) + mass(n)*vxh(:,n)
      enddo

      vxb(:,1) = -vxtmp(:)/msys

      do n=2,nbod
        vxb(:,n) = vxh(:,n) + vxb(:,1)
      enddo

      return
      end     ! coord_vh2b_symbap
c--------------------------------------------------------------------------

