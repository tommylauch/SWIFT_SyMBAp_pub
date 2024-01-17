!*************************************************************************
!                            UTIL_EXIT.F
!*************************************************************************
! Exits program
!             Input:
!                 iflg          ==>  status of exit
!                                       = 0 if normal exit
!                                       = 1 if exit because error
! Remarks: 
! Authors:  Hal Levison 
! Date:    8/6/93
! Last revision: MD : change calc. of rhil Apr. 25

subroutine util_exit(iflg)
use swift_mod
implicit none
integer, intent(in) :: iflg

!...  Executable code 

   write(*,*) ' '
   if(iflg.eq.0) then
      write(*,*) 'Normal termination of SWIFT.'
   else
      write(*,*) 'Terminating SWIFT due to ERROR!!! '
   endif
   write(*,*) '----------------------------------------------------'

stop
end subroutine util_exit
