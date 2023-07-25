c*************************************************************************
c                            SYMBA5P_STEP_INTERP.F
c*************************************************************************
c
c             Input:
c                 time          ==> Current time (real scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 rhill         ==>  Radius of hill sphere (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 lclose        ==> .true. --> marge particles if they
c                                    get too close. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 eoff          ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c                mtiny          ==>  Small mass  (real array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 rpl           ==>  Recalculated physical size of a planet.
c                                    if merger happened (real array)
c                 nbod          ==>  Recalculated number of massive bodies 
c                                    if merger happened (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass          ==>  Recalculated mass of bodies 
c                                    if merger happened (real array)
c                 mergelst      ==>  list of mergers (int array)
c                 mergecnt      ==>  count of mergers (int array)
c                 eoff          ==>  Energy offset (real scalar)
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/21/96
c Last revision: 5/13/99

      subroutine symba5p_step_interp(time,iecnt,ielev,nbod,
     &     nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,
     &     vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,mtiny,
     &     grpie,grppc,grpc)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      real*8 mass(nbod),dt,j2rp2,j4rp4,time,mtiny
      integer iecnt(nbod),ielev(nbod)
      logical*2 lclose 
      integer ielst(2,NENMAX),ielc

c...  Inputs and Outputs:
      integer nbod,nbodm
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 rpl(nbod),eoff,rhill(nbod)

c...  Outputs
      integer mergelst(2,nbod),mergecnt

c...  Internals:
      integer irec,ilevl(NTPMAX),i
      real*8 dth
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      real*8 ptxe,ptye,ptze
      logical*1 svdotr(NENMAX)  ! Used by symba_step_recur

c... grouping
      integer j,k,ip2
      integer gplst(2,GRPMAX)
      integer grpie(GRPMAX,GRPNMAX),grppc(GRPNMAX),grpc
      
c... grouping ^^^^^^^^

      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  Convert vel to bery to jacobi coords
      call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxb,ptyb,ptzb)

c...  Get the accelerations in helio frame. For each object
c...     only include those guys that it is not encountering with. 
      call symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Do a recursion step for full dt for particles not in close encounter
      irec = -1
      call symba5p_helio_drift(nbod,ielev,irec,mass,xh,
     &     yh,zh,vxb,vyb,vzb,dt)
      irec = 0
      do i=2,nbod
         ilevl(i) = 0
      enddo
      mergecnt = 0
c------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(i,j,svdotr,gplst)
!$OMP& SHARED(ielst,time,nbod,nbodm,mass,irec,eoff,
!$OMP& ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,rpl,mergelst,
!$OMP& mergecnt,dt,grpc,grppc,grpie,iecnt,ilevl)
!$OMP DO SCHEDULE(GUIDED)
c...  Guided scheduling specified as the groups have different sizes
c...   and goes to different levels
      do j=1,grpc
c...  cheat symba5_step_recur with gplst and grppc(j)
         do i=1,grppc(j)
             gplst(1,i) = ielst(1,grpie(i,j))
             gplst(2,i) = ielst(2,grpie(i,j))
         enddo
      call symba5p_step_recur(time,nbod,nbodm,mass,irec,ilevl,
     &     iecnt,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,
     &     rpl,mergelst,mergecnt,dt,eoff,svdotr,grppc(j),gplst)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c..   Remove any encounters with ip2s in the encounter matrix
      do i=1,mergecnt
          ip2=mergelst(2,i)
          j = 1
          do while(j.le.ielc)
            if( (ielst(1,j).eq.ip2) .or. (ielst(2,j).eq.ip2) ) then
               do k=j+1,ielc
                  ielst(1,k-1) = ielst(1,k)
                  ielst(2,k-1) = ielst(2,k)
               enddo
               ielc = ielc - 1
            else
               j = j + 1
            endif
          enddo
      enddo
c------------------------------

c...  Get the accelerations in helio frame. For each object
c...     only include those guys that it is not encountering with. 
      call symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxe,ptye,ptze)

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      return
      end   ! symba5p_step_interp
c---------------------------------------------------------------------

