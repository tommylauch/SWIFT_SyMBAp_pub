c*************************************************************************
c                            IO_WRITE_HDR
c*************************************************************************
c write out header part of the real*8 binary file
c
c             Input:
c                 iu              ==> unit number to write to
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 istat           ==>  status of the test paricles
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

subroutine io_write_hdr(iu,time,nbod,ntp,istat) 
implicit none
use swift_mod
use io_mod

integer(ik), intent(in) :: nbod,ntp,istat(NTPMAX,NSTAT),iu
real(rk), intent(in)    :: time

integer(ik)             :: i,nleft

c...  Executable code
c...  calculate number of remaining test particles
nleft = 0
do i=1,ntp
   if (istat(i,1).eq.0) then
      nleft = nleft + 1
   endif
enddo

write(iu) time,nbod,nleft

return
end subroutine io_write_hdr
