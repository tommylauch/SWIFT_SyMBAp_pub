c*************************************************************************
c                        KICKVH_SYMBAP.F
c*************************************************************************
c To kick the velocity components vxh(*) by axh(*)*dt 
c
c             Input:
c                 nbod          ==>  number of bodies (int scalar)
c                 vxh           ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 axh           ==>  acceleration in helio coord
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 vxh           ==>  final velocity in helio coord 
c                                    (real arrays)
c
c     ALGORITHM: Obvious  
*     REMARKS:  Only alters particles 2 thru nbod since Sun is #1
c       
c     AUTHOR:  M. Duncan.
c     DATE WRITTEN:  Feb. 2, 1993.
c     REVISIONS: 2/18/93   HFL


      subroutine kickvh_symbap(nbod,vxh,axh,dt) 


      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 axh(3,nbod),dt

c...   Inputs and Output:
      real*8 vxh(3,nbod)

c...  Internals:
      integer n

c----
c...  Executable code 
      do n=2,nbod
         vxh(:,n) = vxh(:,n) + axh(:,n)*dt
      enddo
      return
      end    ! kickvh_symbap
c-----------------------------------------------------------------------------