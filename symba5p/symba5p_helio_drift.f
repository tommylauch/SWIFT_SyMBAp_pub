c*************************************************************************
c                        SYMBA5P_HELIO_DRIFT.F
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

      subroutine symba5p_helio_drift(nbod,ielev,irec,mass,xh,yh,zh,
     &     vxb,vyb,vzb,dt)	

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer nbod,irec
      real*8 mass(NTPMAX),dt
      integer ielev(NTPMAX)

c...  Inputs and Outputs:
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(j,iflg)
!$OMP& SHARED(nbod,mass,xh,yh,zh,vxb,vyb,vzb,
!$OMP& dt,ielev,irec)
!$OMP DO
      do j = 2,nbod
         if( (ielev(j).eq.irec) .and. (mass(j).ne.0.0d0) ) then
            call drift_one(mass(1),xh(j),yh(j),zh(j),
     &           vxb(j),vyb(j),vzb(j),dt,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
               write(*,*) mass(1),dt
               write(*,*) xh(j),yh(j),zh(j),' H '
               write(*,*) vxb(j),vyb(j),vzb(j),' B '
               write(*,*) ' STOPPING '
               call util_exit(1)
            endif
         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end
c--------------------------------------------------------------------------
