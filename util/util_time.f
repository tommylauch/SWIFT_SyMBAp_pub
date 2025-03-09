!*************************************************************************
!                            UTIL_time.F
!*************************************************************************
! print system time
!             Output:
!                 yr         ==> year
!                 mo         ==> month
!                 day        ==> day
!                 hr         ==> hour
!                 mm         ==> minute
!                 sec        ==> second
!             Internal:
!              values        ==> array of integers for intrinsic subroutine
!                                date_and_time
! Remarks: 
! Authors:  T!H Lau
! Date:    3/8/2024
! Last revision: 

      subroutine util_time(yr,mo,day,hr,mm,sec)
      
      include '../swift.inc'

!...  Outputs:
      integer yr,mo,day,hr,mm,sec

!...  Internals
      integer values(8)
!-----
!...  Executable code 

      call date_and_time(values=values)
      yr = values(1)
      mo = values(2)
      day = values(3)
      hr = values(5)
      mm = values(6)
      sec = values(7)
      return
      end  ! util_time

!---------------------------------------------------
