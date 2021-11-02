c*************************************************************************
c                            IO_READ_MASS
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c read in the mass file.
c
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 iu              ==> unit number to read to
c
c             Returns:
c               io_read_mass     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c
c Remarks: Based on io_read_frame
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      integer function io_read_mass(time,nbod,mass,iu)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Outputs
      integer nbod
      real*8 mass(nbod),time

c...  Internals
      real*4 mass4(NTPMAX)
      real*4 ttmp
      integer i,ierr,ixdrreal,ixdrint
      integer ixdrrmat

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

      ierr = ixdrreal( iu, ttmp )
      io_read_mass = ierr
      if(ierr.ne.0) then
         return
      endif
      time = ttmp

      ierr = ixdrint( iu, nbod )
      io_read_mass = ierr
      if(ierr.ne.0) then
         return
      endif

      ierr = ixdrrmat ( iu, nbod, mass4 )
      io_read_mass = ierr
      if(ierr.ne.0) then
         return
      endif

      do i=1,nbod
         mass(i) = mass4(i)
      enddo

#else

      write(*,*) ' SWIFT ERROR: in io_write_frame: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end      ! io_read_mass
c----------------------------------------------------------------------

