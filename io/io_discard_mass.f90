c*************************************************************************
c                            IO_DISCARD_MASS
c*************************************************************************
c Write out information about a discarded massive body.
c
c             Input:
c                 init          ==>  initiize flag if = 0 initialize and return
c                                                     = 1 run through 
c                 id            ==> particle number (int scalar)
c                 time          ==>  current time (real scalar)
c                 m1            ==>  Mass of pl (real scalar)
c                 r1            ==>  Radius of pl 2 (real scalar)
c                 x1            ==>  current position of pl 1 in helio coord 
c                                    (real scalar)
c                 vx1           ==>  current velocity of pl 1 in helio coord 
c                                    (real scalar)
c                 iu            ==> IO unit (int scalar)
c                 iwhy          ==> reason for discard (int scalar)
c                 fopenstat     ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 9/11/09

subroutine io_discard_mass(init,time,id,m1,r1,x1,vx1,iu,iwhy,fopenstat)
implicit none
use swift_mod
use io_mod

integer(ik), intent(in)        :: iwhy,iu,init,id
real(rk), intent(in)           :: time,m1,r1
real(rk), intent(in)           :: x1(:),vx1(:)
character(len = :), intent(in) :: fopenstat

integer(ik)                    :: ierr
c...  Executable code 
if (init.eq.0) then
   call io_open(iu,'discard_mass.out',fopenstat,'FORMATTED',ierr)
c... if there was an error and fopenstat='append', try to open as new
   if (ierr.ne.0) then  
      if ( (fopenstat(1:6).eq.'append') .or.                           &
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

write(iu,'1x,1p1e23.16,1x,i4') time,iwhy
write(iu,'i7,1x,i7,1x,2(1p1e23.16,1x)') -1,id,m1,r1
write(iu,'3(1p1e23.16,1x)') x1(1),x1(2),x1(3)
write(iu,'3(1p1e23.16,1x)') vx1(1),vx1(2),vx1(3)

close(unit = iu)
return
end subroutine io_discard_mass
