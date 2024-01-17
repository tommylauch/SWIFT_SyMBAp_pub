!*************************************************************************
!                            IO_DISCARD_MASS
!*************************************************************************
! Write out information about a discarded massive body.
!             Input:
!                 init          ==>  initiize flag if = 0 initialize and return
!                                                     = 1 run through 
!                 id            ==> particle number (int scalar)
!                 time          ==>  current time (real scalar)
!                 m1            ==>  Mass of pl (real scalar)
!                 r1            ==>  Radius of pl 2 (real scalar)
!                 x1            ==>  current position of pl 1 in helio coord 
!                                    (real scalar)
!                 vx1           ==>  current velocity of pl 1 in helio coord 
!                                    (real scalar)
!                 iu            ==> IO unit (int scalar)
!                 iwhy          ==> reason for discard (int scalar)
!                 fopenstat     ==>  The status flag for the open 
!                                      statements of the output files.  
!                                          (character*80)
! Remarks: 
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision: 9/11/09

subroutine io_discard_mass(init,time,id,m1,r1,x1,vx1,iu,iwhy,fopenstat)
use swift_mod
use util_interface
use io_interface, except_this_one => io_discard_mass
implicit none

integer(ik), intent(in)        :: iwhy,iu,init,id
real(rk), intent(in)           :: time,m1,r1
real(rk), intent(in)           :: x1(:),vx1(:)
character(len=50), intent(in)  :: fopenstat

integer(ik)                    :: ierr
!...  Executable code 
   if (init.eq.0) then
      call io_open(iu,'discard_mass.out',fopenstat,'FORMATTED',ierr)
!... if there was an error and fopenstat='append', try to open as new
      if (ierr.ne.0) then  
         if ( (fopenstat(1:6).eq.'append') .or.                        &
              (fopenstat(1:6).eq.'APPEND') ) then
            call io_open(iu,'discard_mass.out','new','FORMATTED',ierr)
         endif
      endif

      if (ierr.ne.0) then
         write(*,*) ' SWIFT ERROR: in io_discard_mass: '
         write(*,*) '    Could not open discard output file'
         call util_exit(1)
      endif
      close(unit = iu)
      return
   else
      call io_open(iu,'discard_mass.out','append','FORMATTED',ierr)
   endif

   write(iu,fmt_ri)   time,iwhy
   write(iu,fmt_iirr) -1_ik,id,m1,r1
   write(iu,fmt_rrr)  x1(1),x1(2),x1(3)
   write(iu,fmt_rrr)  vx1(1),vx1(2),vx1(3)

   close(unit = iu)

return
end subroutine io_discard_mass
