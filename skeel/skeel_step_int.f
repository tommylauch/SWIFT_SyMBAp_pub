c*************************************************************************
c                             SKEEL_STEP_INT.F
c*************************************************************************
c Does a step for interacting particles
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of tesp particles (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c              xbeg,ybeg,zbeg   ==>  initial position of planets
c                                    (real arrays)
c           vxbeg,vybeg,vzbeg   ==>  initial position of planets
c                                    (real arrays)
c              xend,yend,zend   ==>  final position of planets
c                                    (real arrays)
c           vxend,vyend,vzend   ==>  final position of planets
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real arrays)
c                 ptxb,ptyb,ptzb ==> beg momentum of sun
c                                       (real scalars)
c                 ptxe,ptye,ptze ==> end momentum of sun
c                                       (real scalars)
c                vxsb,vxsb,vxsb  ==> Initial vel of the Sun: tp's need this
c                                         (real scalars)
c                vxse,vxse,vxse  ==> final vel of the Sun: tp's need this
c                                       (real scalars)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 ienc          ==> ienc(j) = 0 if tp j not involved in enc 
c                                           = planet# if it is. 
c                                           (integer array)
c                 dt            ==>  time step
c                 time          ==>  time
c
c             Output:
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c                                    NOTE: only the tp involed in an encounter
c                                          will have their x and v's changed 
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                    peri=QFLG do not spill particle
c                                         (real array)
c
c
c
c Remarks: Now based on helio
c Authors:  Hal Levison 
c Date:   9/23/96
c Last revision: 2/13/01

      subroutine skeel_step_int(nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,
     &     ptxb,ptyb,ptzb,ptxe,ptye,
     &     ptze,vxsb,vysb,vzsb,vxse,vyse,vzse,     
     &     xht,yht,zht,vxht,vyht,vzht,istat,ienc,dt,time,
     &     isperi,peri)

      include '../swift.inc'
      include 'skeel.inc'

c...  Inputs Only: 
      integer nbod,ntp,ienc(NTPMAX)
      real*8 mass(nbod),dt,time
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 vxbeg(nbod),vybeg(nbod),vzbeg(nbod)
      real*8 xend(nbod),yend(nbod),zend(nbod)
      real*8 vxend(nbod),vyend(nbod),vzend(nbod)
      real*8 j2rp2,j4rp4
      real*8 ptxb,ptyb,ptzb
      real*8 ptxe,ptye,ptze
      real*8 vxsb,vysb,vzsb,vxse,vyse,vzse


c...  Inputs & Outputs Only: 
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Outputs Only
      real*8 peri(NTPMAX)   
      integer isperi(NTPMAX)

c...  Internals: 
      integer i,ipl,irec,irecl,istattmp(NTPMAX),j
      real*8 masst(NPLMAX),dth
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)
      real*8 r2hill(NPLMAX),tperi(NTPMAX)   
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)
      real*8 xrel,yrel,zrel,vxrel,vyrel,vzrel
      logical*2 lperi

      integer i1st              ! =0 first time through; =1 after
      data i1st/0/
      save i1st,r2hill
c----
c...  Executable code 

c...  if first time through, calc hill's shere for the planets
      if(i1st.eq.0) then
         call util_hills(nbod,mass,xbeg,ybeg,zbeg,vxbeg,vybeg,
     &        vzbeg,r2hill)
         i1st = 1
      endif

      dth = 0.5d0*dt

      call coord_vh2b_tp(ntp,vxht,vyht,vzht,vxsb,vysb,vzsb,
     &     vxbt,vybt,vzbt)

      do i=1,ntp
         if( (istat(i,1).eq.0) .and. (ienc(i).ne.0) ) then

            ipl = ienc(i)

c...        apply a linear drift
            do j=1,ntp
               istattmp(j) = 1
            enddo
            istattmp(i) = 0      ! only do the particle in question

            call helio_lindrift_tp(ntp,ptxb,ptyb,ptzb,
     &           dth,istattmp,xht,yht,zht)

c....       apply a kick with all planets but the one involved in encounter
            do j=1,nbod
               masst(j) = mass(j)
            enddo
            masst(ipl) = 0.0d0

c...        apply a kick
            call helio_getacch_tp(nbod,ntp,masst,j2rp2,j4rp4,
     &           xbeg,ybeg,zbeg,xht,yht,zht,istattmp,axht,ayht,azht)
            call kickvh_tp(ntp,vxbt,vybt,vzbt,axht,ayht,azht,
     &           istattmp,dth) 

c....       do the encounter part
            irec = 0
            irecl = 0

c...        Set up stuff to check for peri crossing
            xrel = xht(i) - xbeg(ipl)
            yrel = yht(i) - ybeg(ipl)
            zrel = zht(i) - zbeg(ipl)
            vxrel = vxbt(i) - vxbeg(ipl)
            vyrel = vybt(i) - vybeg(ipl)
            vzrel = vzbt(i) - vzbeg(ipl)
            call util_peri1(0,xrel,yrel,zrel,vxrel,vyrel,
     &          vzrel,mass(ipl),isperi(i),peri(i),tperi(i))

            call skeel_step_recur(mass(1),mass(ipl),r2hill(ipl),
     &         xbeg(ipl),ybeg(ipl),zbeg(ipl),vxbeg(ipl),vybeg(ipl),
     &         vzbeg(ipl),xend(ipl),yend(ipl),zend(ipl),vxend(ipl),
     &         vyend(ipl),vzend(ipl),xht(i),yht(i),zht(i),vxbt(i),
     &         vybt(i),vzbt(i),istat(i,1),istat(i,2),irec,irecl,dt,
     &         isperi(i),peri(i))

c...        apply a kick
            call helio_getacch_tp(nbod,ntp,masst,j2rp2,j4rp4,
     &           xend,yend,zend,xht,yht,zht,istattmp,axht,ayht,azht)
            call kickvh_tp(ntp,vxbt,vybt,vzbt,axht,ayht,azht,
     &           istattmp,dth) 

c...        apply a linear drift
            call helio_lindrift_tp(ntp,ptxe,ptye,ptze,
     &           dth,istattmp,xht,yht,zht)

c...        Put back to helio
            vxht(i) = vxbt(i) - vxse
            vyht(i) = vybt(i) - vyse
            vzht(i) = vzbt(i) - vzse

         endif

      enddo

      return
      end      ! Skeel_step_int.f
c---------------------------------------------------------------------
