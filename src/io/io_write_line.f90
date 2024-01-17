!*************************************************************************
!                            IO_WRITE_LINE
!*************************************************************************
! write out one line to real*8 binary file.
!      Input:
!            iu       ==> unit number to write to
!            a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!        omega    ==> argument of perihelion (real scalar)
!        capm     ==> mean anomoly(real scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/22/94
! Last revision: 

subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm)
use swift_mod
implicit none

integer(ik), intent(in) :: iu,id
real(rk), intent(in)    :: a,e,inc,capom,omega,capm

   write(iu) id,a,e,inc,capom,omega,capm

return
end subroutine io_write_line
