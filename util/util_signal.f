!*************************************************************************
!                            util_signal.F
!*************************************************************************
! This module contains system signal handler
!
      module util_signal
      integer sig_recv
      
      contains
         subroutine util_signal_handler
         implicit none
         write(*,*) 'SIGTERM Received'
         sig_recv = 1
         end subroutine util_signal_handler
      end module util_signal
!---------------------------------------------------------------------