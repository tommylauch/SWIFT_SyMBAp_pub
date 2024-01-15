c*************************************************************************
c                            IO_WRITE_MASS_R
c*************************************************************************
c write out masses
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    1/9/97
c Last revision: 11/2/99

subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)
implicit none
use swift_mod
use io_mod

integer(ik), intent(in)        :: nbod,iu
real(rk), intent(in)           :: mass(:),time
character(len = :), intent(in) :: oname,fopenstat

c...  Internals
integer(ik)                    :: ierr,i
integer(ik), save              :: ldir,lfile
character(len = :), save       :: dirname,filename

integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after

c...  Executable code 
c...  if first time through open file
if(i1st.eq.0) then
   call io_splitname(oname,dirname,ldir,filename,lfile)
   call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),        &
                fopenstat,'UNFORMATTED',ierr)
   if(ierr.ne.0) then
      write(*,*) ' SWIFT ERROR: in io_write_mass_r: '
      write(*,*) '     Could not open binary output file:'
      call util_exit(1)
   endif
   i1st = 1_ik
else
   call io_open(iu,dirname(1:ldir)//'mass.'//filename(1:lfile),        &
                'append','UNFORMATTED',ierr)
endif

write(iu) time,nbod
write(iu) (mass(i),i=1,nbod)

close(iu)
return
end subroutine io_write_mass
