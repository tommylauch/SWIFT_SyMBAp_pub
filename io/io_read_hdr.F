c*************************************************************************
c                            IO_READ_HDR
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c
c             Input:
c                 iu            ==> unit number to write to
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nleft         ==>  number of active tp (int scalar)
c
c             Returns:
c               io_read_hdr     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      integer function io_read_hdr(iu,time,nbod,nleft) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output
      integer nbod,nleft
      real*8 time

c...  Internals
      real*4 ttmp
      integer nn(2)
      integer ierr,ixdrreal,ixdrimat

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

      ierr = ixdrreal( iu, ttmp )
      io_read_hdr = ierr
      if(ierr.ne.0) then
         return
      endif

      ierr = ixdrimat( iu, 2, nn )
      io_read_hdr = ierr
      if(ierr.ne.0) then
         return
      endif

      nbod = nn(1)
      nleft = nn(2)
      time = ttmp

#else

      write(*,*) ' SWIFT ERROR: in io_read_hdr: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end     ! io_read_hdr.f
c---------------------------------------------------------------------------

