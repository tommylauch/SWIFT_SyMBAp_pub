!*************************************************************************
!                            IO_READ_HDR
!*************************************************************************
! read in header part of the real*8 file
!             Input:
!                 iu            ==> unit number to write to
!             Output:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nleft         ==>  number of active tp (int scalar)
!             Returns:
!               io_read_hdr     ==>   =0 read ok
!                                    !=0 read failed is set to iostat variable
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/22/94
! Last revision: 

function io_read_hdr(iu,time,nbod,nleft)
use swift_mod
implicit none
integer(ik)              :: io_read_hdr

integer(ik), intent(in)  :: iu

integer(ik), intent(out) :: nbod,nleft
real(rk), intent(out)    :: time

integer(ik)              :: ierr

   read(iu,iostat=ierr) time,nbod,nleft
   io_read_hdr = ierr
   if (ierr.ne.0_ik) return

return
end function io_read_hdr
