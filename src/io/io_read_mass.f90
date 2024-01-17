!*************************************************************************
!                            IO_READ_MASS_R
!*************************************************************************
! read in the mass file.
!             Output:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 iu              ==> unit number to read to
!             Returns:
!               io_read_mass     ==>   =0 read ok
!                                    !=0 read failed is set to iostat variable
! Remarks: Based on io_read_frame
! Authors:  Hal Levison 
! Date:    1/9/97
! Last revision: 11/2/99

function io_read_mass(time,nbod,mass,iu)
use swift_mod
implicit none
integer(ik)                    :: io_read_mass

integer(ik), intent(in)        :: iu

integer(ik), intent(out)       :: nbod
real(rk), intent(out)          :: mass(:),time

integer(ik)                    :: i,ierr

   read(iu,iostat=ierr) time,nbod
   io_read_mass = ierr
   if (ierr.ne.0_ik) return

   read(iu,iostat=ierr) (mass(i),i=1,nbod)
   io_read_mass = ierr
   if (ierr.ne.0_ik) return

return
end function io_read_mass
