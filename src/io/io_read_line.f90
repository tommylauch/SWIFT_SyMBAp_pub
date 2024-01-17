!*************************************************************************
!                            IO_READ_LINE_R
!*************************************************************************
! read one line from real*4 binary file.
!      Input:
!            iu       ==> unit number to write to
!      Output:
!             a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!        omega    ==> argument of perihelion (real scalar)
!        capm     ==> mean anomoly(real scalar)
!       Returns:
!      io_read_line    ==>   =0 read ok
!                           !=0 read failed is set to iostat variable
! Remarks: 
! Authors:  Hal Levison 
! Date:    2/22/94
! Last revision: 

function io_read_line(iu,id,a,e,inc,capom,omega,capm) 
use swift_mod
implicit none
integer(ik)                    :: io_read_line

integer(ik), intent(in)        :: iu

integer(ik), intent(out)       :: id
real(rk), intent(out)          :: a,e,inc,capom,omega,capm

integer(ik)                    :: ierr

   read(iu,iostat=ierr) id,a,e,inc,capom,omega,capm
   io_read_line = ierr
   if (ierr.ne.0_ik) return

return
end function io_read_line

