!*************************************************************************
!                        DRIFT_KEPU_STUMPFF.F
!*************************************************************************
! subroutine for the calculation of stumpff functions
! see Danby p.172  equations 6.9.15
!             Input:
!                 x             ==>  argument
!             Output:
!                 c0,c1,c2,c3   ==>  c's from p171-172
!                                       (real scalors)
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 2/3/93

subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)
use swift_mod
implicit none

real(rk), intent(inout) :: x

real(rk), intent(out)   :: c0,c1,c2,c3

integer(ik)             :: n,i
real(rk)                :: xm

!...  Executable code 

   n = 0_ik
   xm = 0.1_rk
   do while (abs(x).ge.xm)
      n = n + 1
      x = x*0.25
   enddo

   c2 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/182.)                      &
           /132.)/90.)/56.)/30.)/12.)/2.
   c3 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/210.)                      &
           /156.)/110.)/72.)/42.)/20.)/6.
   c1 = 1. - x*c3
   c0 = 1. - x*c2

   if (n.ne.0) then
      do i=n,1,-1
         c3 = (c2 + c0*c3)/4.
         c2 = c1*c1/2.
         c1 = c0*c1
         c0 = 2.*c0**2 - 1.
         x = x * 4.
       enddo
    endif

return
end subroutine drift_kepu_stumpff
