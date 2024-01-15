c*************************************************************************
c                        HELIO_DRIFT.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
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
c Remarks:  Based on drift.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 1/8/97  for symba

      subroutine helio_drift_symbap(nbod,mass,xh,vxb,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(:),dt

c...  Inputs and Outputs:
      real*8 xh(:,:)
      real*8 vxb(:,:)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(j,iflg)
!$OMP& SHARED(nbod,mass,xh,vxb,dt)
!$OMP DO
      do j=2,nbod
         if(mass(j).ne.0.0d0) then
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
