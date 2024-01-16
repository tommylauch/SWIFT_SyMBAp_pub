!************************************************************************
!                          IO_DUMP_PARAM.F
!************************************************************************
! IO_DUMP_PARAM dumps out the parameters for the integration. 
!      Input:
!       dparfile      ==>  Name of file to write to (character*80)
!            t0       ==> Initial time (real scalar)
!            tstop    ==> final time (real scalar)
!            dt       ==> time step  (real scalar)
!            dtout    ==> time between binary outputs (real scalar)
!            dtdump   ==> time between dumps  (real scalar)
!            iflgchk  ==>  =0 don't run diagnostic routines
!                         !=0 run them
!      rmin,rmax      ==>  maximum and min distance from Sun
!                                if <0  then don't check
!                                    (real scalar)
!       rmaxu         ==>  maximum distance from Sun in not bound
!                                 if <0  then don't check
!                                      (real scalar)
!       qmin          ==> Smallest perihelion distance
!                                 if <0  then don't check
!                                      (real scalar)
!       lclose        ==> .true. --> discard particle if it gets 
!                                    too close to a planet. Read in that 
!                                    distance in io_init_pl
!                                      (logical*2 scalar)
!       outfile       ==>  Name of binary output file (character*80)
! Remarks: 
! Authors:  Martin Duncan
! Date:    3/2/93 
! Last revision:  5/10/94 HFL

subroutine io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,             &
                         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)	
implicit none
use swift_mod
use io_mod
use io_interface, except_this_one => io_dump_param

real(rk), intent(in)           :: t,tstop,dt
integer(ik), intent(in)        :: iflgchk
real(rk), intent(in)           :: dtout,dtdump
real(rk), intent(in)           :: rmin,rmax,rmaxu,qmin
logical(ik), intent(in)        :: lclose
character(len = :), intent(in) :: outfile,dparfile

character(len = :)             :: lflg(0:IO_NBITS-1),cclose
integer(ik)                    :: i,ierr

!...  Executable code 

! Open parameter data file for the dump
   call io_open(7,dparfile,'unknown','formatted',ierr)

	write(7,*) t,tstop,dt
	write(7,*) dtout,dtdump

   do i=0,IO_NBITS-1
      if(btest(iflgchk,i)) then 
         lflg(i) = 'T'
      else
         lflg(i) = 'F'
      endif
   enddo

   write(7,'100(a1,1x)') (lflg(i),i=IO_NBITS-1,0,-1)

   if(btest(iflgchk,4)) then ! bit 4 is set
      if (lclose) then
         cclose = 'T'
      else
         cclose = 'F'
      endif
      write(7,*) rmin,rmax,rmaxu,qmin,' ',cclose
   endif

   if(btest(iflgchk,0).or.btest(iflgchk,1)) then 
      write(7,'a') outfile
   endif

   write(7,'a') 'append'

   close(unit = 7)

return
end subroutine io_dump_param
