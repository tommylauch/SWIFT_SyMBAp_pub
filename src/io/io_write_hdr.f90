!*************************************************************************
!                            IO_WRITE_HDR
!*************************************************************************
! write out header part of the real*8 binary file
!             Input:
!                 iu              ==> unit number to write to
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 ntp            ==>  number of massive bodies (int scalar)
!                 istat           ==>  status of the test paricles
! Remarks:
! Authors:  Hal Levison 
! Date:    2/22/94
! Last revision: 

subroutine io_write_hdr(iu,time,nbod,ntp,istat)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod,ntp,istat(:,:),iu
real(rk), intent(in)    :: time

integer(ik)             :: i,nleft

!...  Executable code
!...  calculate number of remaining test particles
   nleft = 0_ik
   do i=1,ntp
      if (istat(i,1).eq.0) then
         nleft = nleft + 1
      endif
   enddo

   write(iu) time,nbod,nleft

return
end subroutine io_write_hdr
