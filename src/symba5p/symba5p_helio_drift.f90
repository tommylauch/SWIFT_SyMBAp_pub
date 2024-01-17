!*************************************************************************
!                        SYMBA5P_HELIO_DRIFT.F
!*************************************************************************
! This subroutine loops thorugh the particles and calls the danby routine
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ielev         ==>  Level of particles (int array)
!                 irec          ==>  current level of the code
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb           ==>  initial position in bary coord 
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xh            ==>  final position in helio coord 
!                                       (real arrays)
!                 vxb           ==>  final position in bary coord 
!                                       (real arrays)
! Remarks:  Based on helio_drift.f
! Authors:  Hal Levison 
! Date:    1/20.97
! Last revision: 

subroutine symba5p_helio_drift(nbod,ielev,irec,mass,xh,vxb,dt)
use swift_mod
use symba5p_mod
use mvs_interface
use util_interface
implicit none

integer(ik), intent(in) :: nbod,irec,ielev(:)
real(rk), intent(in)    :: mass(:),dt

real(rk), intent(inout) :: xh(:,:),vxb(:,:)

integer(ik)             :: j,iflg

!...  Executable code

! Take a drift forward dth
!$OMP PARALLEL DEFAULT (NONE)                 &
!$OMP PRIVATE(j,iflg)                         &
!$OMP SHARED(nbod,mass,xh,vxb,dt,ielev,irec)
!$OMP DO
   do j=2,nbod
      if( (ielev(j).eq.irec) .and. (mass(j).ne.0.0_rk) ) then
         call drift_one(mass(1),xh(1:3,j),vxb(1:3,j),dt,iflg)
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
end subroutine symba5p_helio_drift
