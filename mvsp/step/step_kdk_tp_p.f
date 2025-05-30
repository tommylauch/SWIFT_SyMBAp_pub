c*************************************************************************
c                            STEP_KDK_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xbeg,ybeg,zbeg ==>  massive part position at beginning of dt
c                                       (real arrays)
c                 xend,yend,zend ==>  massive part position at end of dt
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/24/94

      subroutine step_kdk_tp_p(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,xend,xht,vxht,istat,dt)

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(3,NPLMAX),xend(3,NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(3,ntp),vxht(3,ntp)
c...  Internals:
c      integer i
      real*8 dth 
      real*8 axht(3,NTPMAX)

      save axht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp_p(nbod,ntp,mass,j2rp2,j4rp4,xbeg,
     &                  xht,istat,axht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp_p(ntp,vxht,axht,istat,dth)

c...  Take a drift forward full step
      call drift_tp_p(ntp,mass(1),xht,vxht,dt,istat)

c...  Get the accelerations in helio frame.
      call getacch_tp_p(nbod,ntp,mass,j2rp2,j4rp4,xend,
     &                  xht,istat,axht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp_p(ntp,vxht,axht,istat,dth)

      return
      end   ! step_kdk_tp
c---------------------------------------------------------------------

