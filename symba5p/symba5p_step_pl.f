c*************************************************************************
c                            SYMBA5P_STEP_PL.F
c*************************************************************************
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  location of the last massie body 
c                                    (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh            ==>  initial position in helio coord
c                                    (real arrays)
c                 vxh           ==>  initial velocity in helio coord
c                                    (real arrays)
c                 dt            ==>  time step
c                 lclose        ==> .true. --> marge particles if they
c                                    get too close. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 eoff          ==>  Energy offset (real scalar)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c                 mtiny         ==>  Small mass  (real array)
c             Output:
c                 xh            ==>  final position in helio coord
c                                       (real arrays)
c                 vxh           ==>  final velocity in helio coord
c                                       (real arrays)
c                 rpl           ==>  Recalculated physical size of a planet.
c                                    if merger happened (real array)
c                 nbod          ==>  Recalculated number of massive bodies 
c                                    if merger happened (int scalar)
c                 mass          ==>  Recalculated mass of bodies 
c                                    if merger happened (real array)
c                 isenc         ==>  0 --> No encounter during last dt
c                                    1 --> There was encounters
c                                     (integer scalar)
c                 mergelst      ==>  list of mergers (int array)
c                 mergecnt      ==>  count of mergers (int array)
c                 iecnt         ==>  Number of encounters (int*2 array)
c                 eoff          ==>  Energy offset (real scalar)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c                 grpie(i,j)    ==>  lists of pairs in groups (2D integer array)
c                                    stores the index ie in ielst for pair no. i 
c                                    in group no. j
c                 grppc(i)      ==>  Number of pairs in group i (integer array)
c                 grpc          ==>  Number of groups (integer)
c
c Remarks: Based on symba2_step_pl.f 
c Authors:  Hal Levison
c Date:    11/27/97
c Last revision: 

      subroutine symba5p_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,vxh,dt,lclose,rpl,isenc,
     &     mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(*),dt,time,j2rp2,j4rp4,mtiny
      logical*2 lclose

c...  Inputs and Outputs:
      real*8 xh(3,*),vxh(3,*)
      real*8 rpl(*),eoff,rhill(*)

c...  Outputs only
      integer isenc
      integer iecnt(*)
      integer mergelst(2,*),mergecnt

c...  Internals
      integer i,j,ieflg,irec,ij_tm,ij
      logical*1 svdotr            ! Not used in the routine
      integer ielst(2,NENMAX),ielc,ielev(NTPMAX)

c... grouping
      integer grpie(GRPMAX,GRPNMAX),grppc(GRPNMAX),grpc
c----
c...  Executable code 

      do i=2,nbod
         iecnt(i) = 0
         ielev(i) = -1
      enddo
      do i=1,GRPNMAX
         grppc(i) = 0
      enddo
      isenc = 0
c...  check for encounters
      irec = 0
      ielc = 0
      ij_tm = ((nbodm-1)*(nbodm-2))/2
      grpc = 0
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(i,j,ij,ieflg,svdotr)
!$OMP& SHARED(ij_tm,rhill,nbod,nbodm,mass,xh,vxh,
!$OMP& dt,irec,iecnt,ielev,ielc,ielst,isenc,grpc,grppc,grpie)
!$OMP DO COLLAPSE(2)
      do j=2,nbodm
         do i=nbodm+1,nbod
            ieflg = 0
            call symba5_chk_symbap(rhill,nbod,i,j,mass,xh,vxh,
     &           dt,irec,ieflg,svdotr)
            if(ieflg.ne.0) then
!$OMP CRITICAL (ENC)
               isenc = 1
               iecnt(i) = iecnt(i) + 1
               iecnt(j) = iecnt(j) + 1
               ielev(i) = 0
               ielev(j) = 0
               ielc = ielc + 1
               if(ielc.gt.NENMAX) then
                  write(*,*) 'ERROR: encounter matrix is filled.'
                  write(*,*) 'STOPPING'
                  call util_exit(1)
               endif
               ielst(1,ielc) = i
               ielst(2,ielc) = j
               call symba5p_group(ielst,ielc,i,j,grpie,grppc,grpc)
!$OMP END CRITICAL (ENC)
            endif
         enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(GUIDED)
c...  Guided scheduling specified as iterations with close encounter are longer
      do ij=0,(ij_tm-1)
         j=ij/(nbodm-2)+2 !j goes from 2
         i=mod(ij,(nbodm-2))+3  !i goes from 3
         if(i.le.j) then
            j = nbodm-j+2 !j goes to nbodm-1
            i = nbodm-i+3 !i goes to nbodm
         endif
         ieflg = 0
         call symba5_chk_symbap(rhill,nbod,i,j,mass,xh,vxh,
     &           dt,irec,ieflg,svdotr)
         if(ieflg.ne.0) then
!$OMP CRITICAL (ENC)
            isenc = 1
            iecnt(i) = iecnt(i) + 1
            iecnt(j) = iecnt(j) + 1
            ielev(i) = 0
            ielev(j) = 0
            ielc = ielc + 1
            if(ielc.gt.NENMAX) then
               write(*,*) 'ERROR: encounter matrix is filled.'
               write(*,*) 'STOPPING'
               call util_exit(1)
            endif
            ielst(1,ielc) = i
            ielst(2,ielc) = j
            call symba5p_group(ielst,ielc,i,j,grpie,grppc,grpc)
!$OMP END CRITICAL (ENC)
         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL
c...  do a step
      if(isenc.eq.0) then
         call symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,
     &                           xh,vxh,dt)
         mergecnt=0
      else
         call symba5p_step_interp(time,iecnt,ielev,nbod,nbodm,mass,
     &     rhill,j2rp2,j4rp4,lclose,rpl,xh,vxh,dt,mergelst,mergecnt,
     &     eoff,ielc,ielst,mtiny,grpie,grppc,grpc)
         i1st = 0  
      endif

      return
      end ! symba5p_step_pl.f
c-----------------------------------------------------------

