c*************************************************************************
c                            UTIL_EXIT.F
c*************************************************************************
c Exits program
c
c             Input:
c                 iflg          ==>  status of exit
c                                       = 0 if normal exit
c                                       = 1 if exit because error
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    8/6/93
c Last revision: MD : change calc. of rhil Apr. 25

subroutine util_exit(iflg)
implicit none
use swift_mod
integer, intent(in) ::  iflg

c...  Executable code 

write(*,*) ' '

if(iflg.eq.0) then
   write(*,*) 'Normal termination of SWIFT.'
else
   write(*,*) 'Terminating SWIFT due to ERROR!!! '
endif

write(*,*) '----------------------------------------------------'

stop
end subroutine util_exit
