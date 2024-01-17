!*************************************************************************
!                            IO_OPEN.F
!*************************************************************************
! open files
!             Input:
!                 iu              ==>  unit number (integer scalar)
!                 fname           ==>  file name (character*80)
!                 fopenstat       ==>  The status flag for the open 
!                                      statements of the output files.  
!                                          (character*80)
!                 format          ==>  format string (character*80)
!             Output:
!                 ierr            ==>  output from iostat
! Remarks: 
! Authors:  Hal Levison 
! Date:    3/3/94
! Last revision: 1/30/98

subroutine io_open(iu,fname,fopenstat,format,ierr)
use swift_mod
implicit none

integer(ik), intent(in)        :: iu
character(len=50), intent(in)  :: fname,fopenstat,format

integer(ik), intent(out)       :: ierr

!...  Executable code

   if( (fopenstat(1:6).eq.'append') .or. (fopenstat(1:6).eq.'APPEND') )&  
      then
      open(unit=iu,file=fname,status='old',access='append',            &
           form=format,iostat=ierr)
      if(ierr.ne.0) then
         write(*,*) 'Warning:  Could not open ',fname,' with'
         write(*,*) '          position=append.'
         write(*,*) '          Will open as status=new'
         open(unit=iu,file=fname,status='new',form=format,iostat=ierr)
      endif
   else
      open(unit=iu,file=fname,status=fopenstat,form=format,iostat=ierr)
   endif

return
end subroutine io_open
