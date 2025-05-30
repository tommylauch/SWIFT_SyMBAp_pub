c*************************************************************************
c                            STEP_KDK_P.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,           ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh           ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht           ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht           ==>  final position in helio coord 
c                                       (real arrays)
c
c
c Remarks: Adopted from martin's nbwh.f program
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 9/26/97

      subroutine step_kdk_p(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,vxh,xht,vxht,istat,rstat,dt)

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(3,nbod),vxh(3,nbod),xht(3,ntp),vxht(3,ntp)

c...  Internals
      integer i1sttp,i
      real*8 xbeg(3,nbod),xend(3,nbod)

c----
c...  Executable code 

      i1sttp = i1st

c...  remember the current position of the planets
      xbeg = xh

c...  first do the planets
      call step_kdk_pl_p(i1st,nbod,mass,j2rp2,j4rp4,xh,vxh,dt)

      if(ntp.ne.0) then
c...     now remember these positions
         xend = xh

c...     next the test particles
         call step_kdk_tp_p(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &        xbeg,xend,xht,vxht,istat,dt)
      endif

      return
      end   ! step_kdk
c------------------------------------------------------------------------