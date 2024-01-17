!*************************************************************************
!                            IO_WRITE_MASS_R
!*************************************************************************
! write out masses
!             Input:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 oname           ==> output file name (character string) 
!                 iu              ==> unit number to write to
!                 fopenstat       ==>  The status flag for the open 
!                                      statements of the output files.  
! Remarks: Based on io_write_frame
! Authors:  Hal Levison 
! Date:    1/9/97
! Last revision: 11/2/99

subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)
use swift_mod
use util_interface
use io_interface, except_this_one => io_write_mass
implicit none

integer(ik), intent(in)        :: nbod,iu
real(rk), intent(in)           :: mass(:),time
character(len=50), intent(in)  :: oname,fopenstat

integer(ik)                    :: ierr,i
integer(ik), save              :: ldir,lfile
character(len=50), save        :: dirname,filename
integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after

!...  Executable code
!...  if first time through open file
   if(i1st.eq.0) then
      call io_splitname(oname,dirname,ldir,filename,lfile)
      call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),     &
                   fopenstat,'UNFORMATTED',ierr)
      if(ierr.ne.0) then
         write(*,*) ' SWIFT ERROR: in io_write_mass_r: '
         write(*,*) '     Could not open binary output file:'
         call util_exit(1)
      endif
      i1st = 1_ik
   else
      call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),     &
                   'append','UNFORMATTED',ierr)
   endif
   write(iu) time,nbod
   write(iu) (mass(i),i=1,nbod)
   close(iu)

return
end subroutine io_write_mass
