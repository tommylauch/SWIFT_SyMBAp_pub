!*************************************************************************
!                          IO_ENERGY_WRITE.F
!*************************************************************************
! Does the write for anal_jacobi_write
!      Input:
!            i1st           ==>  =0 if first write, =1 if not (int scalar)
!            t              ==>  current time (real scalar)
!            energy         ==>  Total energy
!            eltot          ==>  components of total angular momentum
!                               (real array)
!            iu             ==>  unit to write to
!            fopenstat      ==>  The status flag for the open 
!                                statements of the output files.  
!                                          (character*80)
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/21/94
! Last revision: 3/4/94

subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
use swift_mod
use io_interface, except_this_one => io_energy_write
use util_interface
implicit none

integer(ik), intent(in)        :: iu,i1st
real(rk), intent(in)           :: t,energy,eltot(:)
character(len=50), intent(in)  :: fopenstat

integer(ik)                    :: ierr

!...  Executable code

   if (i1st.eq.0) then
      call io_open(iu,'energy.out',fopenstat,'FORMATTED',ierr)
      if(ierr.ne.0) then
         write(*,*) ' SWIFT ERROR: in anal_energy_write '
         write(*,*) '     Could not open energy.out '
         call util_exit(1)
      endif
   else
      call io_open(iu,'energy.out','append','FORMATTED',ierr)
   endif

   write(iu,'(1x,1p1e12.5,4(2x,1p1e23.16))') t,energy,eltot
   close(iu)

return
end subroutine io_energy_write
