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
c                 xh            ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb           ==>  initial position in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh            ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb           ==>  final position in bary coord 
c                                       (real arrays)
c
c Remarks:  Based on helio_drift.f
c Authors:  Hal Levison 
c Date:    1/20.97
c Last revision: 

      subroutine symba5p_helio_drift(nbod,ielev,irec,mass,xh,vxb,dt)	

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer nbod,irec
      real*8 mass(:),dt
      integer ielev(:)

c...  Inputs and Outputs:
      real*8 xh(:,:),vxb(:,:)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(j,iflg)
!$OMP& SHARED(nbod,mass,xh,vxb,dt,ielev,irec)
!$OMP DO
      do j=2,nbod
         if( (ielev(j).eq.irec) .and. (mass(j).ne.0.0d0) ) then
            call drift_one_symbap(mass(1),xh(1:3,j),vxb(1:3,j),dt,iflg)
            if(iflg.ne.0) then
               write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
               write(*,*) mass(1),dt
               write(*,*) xh(1,j),xh(2,j),xh(3,j),' H '
               write(*,*) vxb(1,j),vxb(2,j),vxb(3,j),' B '
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
