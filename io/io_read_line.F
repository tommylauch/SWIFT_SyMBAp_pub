c*************************************************************************
c                            IO_READ_LINE
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c
c      Input:
c            iu       ==> unit number to write to
c      Output:
c            id       ==> Particle id number (int scalar)
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c       Returns:
c      io_read_line    ==>   =0 read ok
c                           !=0 read failed is set to iostat variable
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      integer function io_read_line(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output: 
      integer id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      real*4 orbel(6)
      integer ierr
      integer ixdrint,ixdrrmat

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

      ierr = ixdrint (iu, id)
      io_read_line = ierr
      if(ierr.ne.0) then
         return
      endif

      ierr = ixdrrmat ( iu, 6, orbel )
      io_read_line = ierr
      if(ierr.ne.0) then
         return
      endif

      a = orbel(1)
      e = orbel(2)
      inc = orbel(3)
      capom = orbel(4)
      omega = orbel(5)
      capm = orbel(6)

#else

      write(*,*) ' SWIFT ERROR: in io_write_line: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end      ! io_read_line
c--------------------------------------------------------------------------

