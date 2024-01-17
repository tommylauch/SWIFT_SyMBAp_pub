!*************************************************************************
!                        SYMBA5P_HELIO_DRIFT_G.F
!*************************************************************************
! This subroutine loops thorugh the particles and calls the danby routine
!             Input:
!                 ielev         ==>  Level of particles (int array)
!                 irec          ==>  current level of the code
!                 mass          ==>  mass of bodies (real array)
!                 xh,yh,zh      ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb,vyb,vzb   ==>  initial position in bary coord 
!                                    (real arrays)
!                 dt            ==>  time step
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
!             Output:
!                 xh,yh,zh      ==>  final position in helio coord 
!                                       (real arrays)
!                 vxb,vyb,vzb   ==>  final position in bary coord 
!                                       (real arrays)
! Remarks:  Based on helio_drift.f
! Authors:  Hal Levison 
! Date:    1/20.97
! Last revision: 

subroutine symba5p_helio_drift_g(ielev,irec,mass,xh,vxb,dt,ielc,ielst)
use swift_mod
use symba5p_mod
use mvs_interface
use util_interface
implicit none

integer(ik), intent(in) :: irec
real(rk), intent(in)    :: mass(:),dt
integer(ik), intent(in) :: ielev(:),ielst(:,:),ielc

real(rk), intent(inout) :: xh(:,:),vxb(:,:)

integer(ik)             :: i,j,iflg,i_ie,j_ie
integer(ik)             :: gpmb(GRPMAX),gpmbc
logical(ik)             :: dup_i,dup_j

!...  Executable code

! Take a drift forward dth
   gpmbc = 0_ik
   do i=1,ielc
      i_ie = ielst(1,i)
      j_ie = ielst(2,i)
      dup_i = .false.
      dup_j = .false.
      do j=1,gpmbc
         if (i_ie.eq.gpmb(j)) then
          dup_i = .true.
         endif
         if (j_ie.eq.gpmb(j)) then
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
      if ( (ielev(j).eq.irec) .and. (mass(j).ne.0.0_rk) ) then
         call drift_one(mass(1),xh(1:3,j),vxb(1:3,j),dt,iflg)
         if(iflg.ne.0) then
            write(*,*) ' Planet ',j,' is lost !!!!!!!!!'
            write(*,*) mass(1),dt
            write(*,*) xh(1,j),xh(2,j),xh(3,j),' H '
            write(*,*) vxb(1,j),vxb(2,j),vxb(3,j),' B '
            write(*,*) ' STOPPING G'
            call util_exit(1)
         endif
      endif
   enddo

return
end subroutine symba5p_helio_drift_g
