c*************************************************************************
c                            IO_WRITE_HDR
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c 
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
c Date:    11/2/99
c Last revision: 

      subroutine io_write_hdr(iu,time,nbod,ntp,istat) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX,NSTAT),iu
      real*8 time

c...  Internals
      integer i,ixdrreal,ixdrimat,ierr
      real*4 ttmp
      integer nleft,nn(2)

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

c...  calculate number of remaining test particles
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
      enddo

      ttmp = time

      ierr = ixdrreal( iu, ttmp )
      if(ierr.gt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_hdr: '
         write(*,*) '     Could not write time'
         call util_exit(1)
      endif
 
      nn(1) = nbod
      nn(2) = nleft
      ierr = ixdrimat( iu, 2, nn )
      if(ierr.gt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_hdr: '
         write(*,*) '     Could not write nbod and nleft'
         call util_exit(1)
      endif

#else

      write(*,*) ' SWIFT ERROR: in io_write_hdr: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end     ! io_write_hdr.f
c---------------------------------------------------------------------------




