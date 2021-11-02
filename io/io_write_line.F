c*************************************************************************
c                            IO_WRITE_LINE
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c
c      Input:
c            iu       ==> unit number to write to
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu,id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      real*4 orbel(6)
      integer ixdrint,ixdrrmat,ierr

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

      ierr = ixdrint (iu, id)
      if(ierr.gt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_line: '
         write(*,*) '     Could not write id'
         call util_exit(1)
      endif

      orbel(1) = a
      orbel(2) = e
      orbel(3) = inc
      orbel(4) = capom
      orbel(5) = omega
      orbel(6) = capm

      ierr = ixdrrmat ( iu, 6, orbel )
      if(ierr.gt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_line: '
         write(*,*) '     Could not write orbit elements',id
         call util_exit(1)
      endif

#else

      write(*,*) ' SWIFT ERROR: in io_write_line: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end      ! io_write_line
c--------------------------------------------------------------------------
