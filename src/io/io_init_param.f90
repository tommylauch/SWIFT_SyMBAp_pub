!************************************************************************
!                          IO_INIT_PARAM.F
!************************************************************************
! INIT_PARAM reads in the parameters for the integration. 
!      Input:
!            infile   ==> File name to read from (character*80)
!      Output:
!            t0       ==> Initial time (real scalar)
!            tstop    ==> final time (real scalar)
!            dt       ==> time step  (real scalar)
!            dtout    ==> time between binary outputs (real scalar)
!            dtdump   ==> time between dumps  (real scalar)
!            iflgchk  ==>  =0 don't run diagnostic routines
!                          bit 0 set ==>  write int*2 binary data file
!                          bit 1 set ==>  write real*4 binary file 
!                          bit 2 set ==>  calc energy of system wrt time
!                          bit 3 set ==>  calc jacobi of the test particles
!                          bit 4 set ==>  check if particles are removed
!                          bit 5 set ==>  include J2 and J4 terms
!      rmin,rmax      ==>  maximum and min distance from Sun
!                                if <0  then don't check
!                                    (real scalar)
!      rmaxu          ==>  maximum distance from Sun in not bound
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
!       fopenstat     ==>  The status flag for the open statements of the
!                          output files.  Must be one of the following:
!                                 new      (die if the file exists)
!                                 append   (add to what is there)
!                                 unknown  (just write over what is there)
! Remarks: 
! Authors:  Martin Duncan
! Date:    3/2/93 
! Last revision:  5/10/94  HFL
subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,iflgchk,      &
                         rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)
use swift_mod
use io_mod
use util_interface
use io_interface, except_this_one => io_init_param
implicit none

character(len=*), intent(in)    :: infile

integer(ik), intent(out)        :: iflgchk
real(rk), intent(out)           :: t0,tstop,dt
real(rk), intent(out)           :: dtout,dtdump
real(rk), intent(out)           :: rmin,rmax,rmaxu,qmin
logical(ik), intent(out)        :: lclose
character(len=*), intent(out)   :: outfile,fopenstat

logical(ik)                     :: lflg(0:IO_NBITS-1)
integer(ik)                     :: i,ierr

!...  Executable code 

   write(*,*) 'Parameter data file is ',infile
! Open and read in parameters
   call io_open(7,infile,'old','formatted',ierr)
   read(7,*) t0,tstop,dt
   write(*,*) 't0,tstop,dt : ',t0,tstop,dt
   read(7,*) dtout,dtdump
   write(*,*) 'dtout,dtdump : ',dtout,dtdump
   read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

   iflgchk=0
   do i=0,IO_NBITS-1
      if(lflg(i)) then
         iflgchk = ibset(iflgchk,i)
      endif
   enddo

   write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

   if (btest(iflgchk,0) .and. btest(iflgchk,1)) then 
      write(*,*) ' SWIFT ERROR: in io_init_param:'
      write(*,*) '    Invalid logical flags '
      write(*,*) '    You cannot request that both a real and '
      write(*,*) '    an integer binary file be written.'
      call util_exit(1)
   endif

   if (btest(iflgchk,4)) then ! bit 4 is set
      read(7,*) rmin,rmax,rmaxu,qmin,lclose
      write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',rmin,rmax,rmaxu,qmin,lclose
   else
      rmin = -1.0_rk
      rmax = -1.0_rk
      rmaxu = -1.0_rk
      qmin = -1.0_rk
      lclose = .false.
   endif

   if (btest(iflgchk,0) .or. btest(iflgchk,1)) then 
      read(7,'(a)') outfile
      write(*,*) 'outfile : ', outfile
      write(*,*) ' '
   endif

   read(7,'(a)') fopenstat
   if ( (fopenstat(1:3).ne.'new') .and.                                &
        (fopenstat(1:3).ne.'NEW') .and.                                &
        (fopenstat(1:7).ne.'unknown') .and.                            &
        (fopenstat(1:7).ne.'UNKNOWN') .and.                            &
        (fopenstat(1:6).ne.'append') .and.                             &
        (fopenstat(1:6).ne.'APPEND') ) then
      write(*,*) ' SWIFT ERROR: in io_init_param:'
      write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
      call util_exit(1)
   endif
        
   close(unit = 7)

return
end subroutine io_init_param
