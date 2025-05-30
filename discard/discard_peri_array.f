c*************************************************************************
c                            DISCARD_PERI_ARRAY.F
c*************************************************************************
c This subroutine checks to see if a partical should be discarded because
c of its perihelion distance gets too small
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 ntp           ==>  number of test bodies (int scalar)
c                 xht,yht,zht    ==>   part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>   part vel in helio coord 
c                                      (real arrays)
c                 qmin            ==>  Smallest perihelion distance 
c                                      (real scalar)
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                 nbod            ==>  Number of planets (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c             Output:
c                 istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 1 if discarded
c	                               istat(i,2) =  -4
c                 rstat           ==>  status of the test paricles
c                                      (2d  real array)
c                                      rstat(i,2) perihelion distance.
c                                      rstat(i,1) time of discard.
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    5/10/94
c Last revision: 1/20/97

         subroutine discard_peri_array(time,ntp,xht,vxht,
     &       qmin,istat,rstat,nbod,mass,xh,vxh)

      include '../swift.inc'

c...  Inputs: 
      integer ntp,nbod
      real*8 mass(nbod),time,qmin
      real*8 xht(3,ntp),vxht(3,ntp)
      real*8 xh(3,nbod),vxh(3,nbod)

c...  Input and Output
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)

c...  internal
      integer i,isperi(NTPMAX),i1st,j
      real*8 peri(NTPMAX),ih,r2
      logical*2 lperi(NTPMAX)
      real*8 r2hill(NTPMAX)

      data i1st/0/
      save i1st,isperi,r2hill

c-----
c...  Executable code 

      if(i1st.eq.0) then     ! if first time through, set things up
         call util_hills_array(nbod,mass,xh,vxh,r2hill)
         call util_peri_array(0,ntp,xht,vxht,mass(1),isperi,peri,lperi)
         i1st = 1
         return                 !  <==== RETURN
      endif

      call util_peri_array(1,ntp,xht,vxht,mass(1),isperi,peri,lperi)

      do i=1,ntp
         if( (istat(i,1).eq.0).and. (isperi(i).eq.0) ) then
            
            ih = 0
            do j=2,nbod
               r2 = sum((xht(:,i)-xh(:,j))**2)
               if(r2.le.r2hill(j)) then
                  ih = 1
               endif
            enddo

            if(ih.eq.0) then
               rstat(i,2) = peri(i)
               if(peri(i).le.qmin) then
                  write(*,*) 'Particle',i,' perihelion distance too',
     &                 ' small at t=',time
                  istat(i,1) =  1
                  istat(i,2) = -4
               endif
            endif

         endif
      enddo

      return
      end       ! discard_peri_array










