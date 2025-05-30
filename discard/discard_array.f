c*************************************************************************
c                            DISCARD_ARRAY.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its position or becuase it becomes unbound
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>   part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>   velocity in helio coord 
c                                        (real arrays)
c                 rmin,rmax      ==>  maximum and min distance from Sun
c                                     if <0  then don't check
c                                        (real scalar)
c                 rmaxu          ==>  maximum distance from Sun in not bound
c                                     if <0  then don't check
c                                        (real scalar)
c                  qmin          ==> Smallest perihelion distance
c                                      if <0  then don't check
c                                          (real scalar)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) closest approach to a planet
c
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c                                      istat(i,2) = -2  -> a<0 & r_helio>rmaxu
c	                               istat(i,2) = -3     r_helio>rmax
c	                               istat(i,2) = -4     q<qmin
c	                               istat(i,2) =  1     r_helio<rmin
c	                               istat(i,2) =  n    too close to planet n
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,2) closest approach to a planet
c                                          as determined by encounter routines.
c
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c	                               istat(i,2) =  n    too close to planet n
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,1) time of discard.
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/2/93
c Last revision: 1/20/97

      subroutine discard_array(time,dt,nbod,ntp,mass,xh,vxh,
     &     xht,vxht,rmin,rmax,rmaxu,qmin,lclose,rplsq,istat,rstat)


      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp
      real*8 mass(nbod),xh(3,nbod),time,rplsq(NPLMAX)
      real*8 vxh(3,nbod),xht(3,ntp),vxht(3,ntp)
      real*8 rmin,rmax,rmaxu,dt,qmin
      logical*2 lclose

c...  Input and Output
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)

c...  Internals
      real*8 xb(3,NPLMAX),vxb(3,NPLMAX)
      real*8 xbt(3,NTPMAX),vxbt(3,NTPMAX),msys

c-----
c...  Executable code 

      if( (rmin.ge.0.0) .or. (rmax.ge.0.0) .or. (rmaxu.ge.0.0) ) then
         call coord_h2b_array(nbod,mass,xh,vxh,xb,vxb,msys)
         call coord_h2b_tp_array(ntp,xht,vxht,xb(:,1),vxb(:,1),xbt,vxbt)
         call discard_sun_array(time,ntp,msys,xht,xbt,
     &       vxbt,rmin,rmax,rmaxu,istat,rstat)
      endif

      if(qmin.ge.0.0) then
         call discard_peri_array(time,ntp,xht,vxht,
     &        qmin,istat,rstat,nbod,mass,xh,vxh)
      endif

      if(lclose) then
         call discard_pl_array(time,dt,nbod,ntp,mass,xh,vxh,
     &       xht,vxht,rplsq,istat,rstat)
      endif

      return
      end    ! discard.f
c-----------------------------------------------------------------------


