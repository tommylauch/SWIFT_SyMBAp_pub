c*************************************************************************
c                        SYMBA5P_HELIO_DRIFT_G.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ielev         ==>  Level of particles (int array)
c                 irec          ==>  current level of the code
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial position in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb,vyb,vzb   ==>  final position in bary coord 
c                                       (real arrays)
c
c Remarks:  Based on helio_drift.f
c Authors:  Hal Levison 
c Date:    1/20.97
c Last revision: 

      subroutine symba5p_helio_drift_g(nbod,ielev,irec,mass,xh,yh,zh,
     &     vxb,vyb,vzb,dt,ielc,ielst)	

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer nbod,irec
      real*8 mass(NTPMAX),dt
      integer ielev(NTPMAX)
      integer ielst(2,NENMAX),ielc
c...  Inputs and Outputs:
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
      
c...  Internals:
      integer i,j,iflg,i_ie,j_ie
      integer gpmb(GRPMAX),gpmbc
      logical*1 dup_i,dup_j
c----
c...  Executable code 

c Take a drift forward dth
      gpmbc = 0
      do i=1,ielc
         i_ie = ielst(1,i)
         j_ie = ielst(2,i)
         dup_i = .false.
         dup_j = .false.
         do j=1,gpmbc
            if(i_ie.eq.gpmb(j))then
             dup_i = .true.
            endif
            if(j_ie.eq.gpmb(j))then
             dup_j = .true.
            endif
         enddo

         if(dup_i.eqv..false.)then
              gpmbc = gpmbc + 1
              gpmb(gpmbc) = i_ie
         endif
         if(dup_j.eqv..false.)then
              gpmbc = gpmbc + 1
              gpmb(gpmbc) = j_ie
         endif
      enddo

      do i=1,gpmbc
         j = gpmb(i)
         if( (ielev(j).eq.irec) .and. (mass(j).ne.0.0d0) ) then
            call drift_one(mass(1),xh(j),yh(j),zh(j),
     &           vxb(j),vyb(j),vzb(j),dt,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
               write(*,*) mass(1),dt
               write(*,*) xh(j),yh(j),zh(j),' H '
               write(*,*) vxb(j),vyb(j),vzb(j),' B '
               write(*,*) ' STOPPING G'
               call util_exit(1)
            endif
         endif
      enddo
      return
      end
c--------------------------------------------------------------------------
