c***********************************************************************
c                          COORD_VB2H_SYMBAP.F
c***********************************************************************
*     PURPOSE: Converts from Barycentric to Helio coords.
*               Velocity only
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*                   mass(*) ==>  masses (real array)
*                                 NOT USED BUT INCLUDED IN ORDER TO HAVE
*                                 SYMMETRY IN SUBROUTINE CALLS
*                 vxb(*) ==> Barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: Obvious 
*     Authors:  Hal Levison
*     WRITTEN:  11/14/96
*     REVISIONS: 11/21/96

      subroutine coord_vb2h_symbap(nbod,mass,vxb,vxh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),vxb(3,nbod)

c...  Outputs:
      real*8 vxh(3,nbod)

c...  Internals:
      integer i

c----
c...  Executable code 

      vxb(:,1) = -mass(2)*vxb(:,2)

      do i=3,nbod
         vxb(:,1) = vxb(:,1) - mass(i)*vxb(:,i)
      enddo

      vxb(:,1) = vxb(:,1)/mass(1)

      do i=2,nbod
         vxh(:,i) = vxb(:,i) - vxb(:,1)
      enddo

      return
      end     ! coord_vb2h_symbap

c--------------------------------------------------------------------------

