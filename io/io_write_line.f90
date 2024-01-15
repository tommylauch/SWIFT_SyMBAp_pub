c*************************************************************************
c                            IO_WRITE_LINE
c*************************************************************************
c write out one line to real*4 binary file.
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
c Date:    2/22/94
c Last revision: 

subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm) 
implicit none
use swift_mod
use io_mod

c...  Inputs: 
integer(ik), intent(in) :: iu,id
real(rk), intent(in)    :: a,e,inc,capom,omega,capm

write(iu) id,a,e,inc,capom,omega,capm

return
end subroutine io_write_line
