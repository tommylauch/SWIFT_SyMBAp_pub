c*************************************************************************
c                            IO_WRITE_FRAME
c*************************************************************************
c write out a whole frame to an real*4 binary file.
c both massive and test particles
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  current position in helio coord 
c                                    (real array)
c                 vxh           ==>  current velocity in helio coord 
c                                    (real array)
c                 xht           ==>  current part position in helio coord 
c                                    (real array)
c                 vxht          ==>  current velocity in helio coord 
c                                    (real array)
c                 istat         ==>  status of the test paricles
c                                    (2d integer array)
c                                    istat(i,1) = 0 ==> active:  = 1 not
c                                    istat(i,2) = -1 ==> Danby did not work
c                 oname         ==>  output file name (character string) 
c                 iu            ==>  unit number to write to
c                 fopenstat     ==>  The status flag for the open 
c                                    statements of the output files.  
c                                    (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

subroutine io_write_frame(time,nbod,ntp,mass,xh,vxh,xht,vxht,istat,    &
                          oname,iu,fopenstat)
implicit none
use swift_mod
use io_mod

integer(ik), intent(in)        :: nbod,ntp,iu
real(rk), intent(in)           :: mass(nbod),time
integer(ik), intent(in)        :: istat(NTPMAX,NSTAT)
real(rk), intent(in)           :: xh(:,:),vxh(:,:)
real(rk), intent(in)           :: xht(:,:),vxht(:,:)
character(len = :), intent(in) :: oname,fopenstat

integer(ik)                    :: i,id
integer(ik)                    :: ialpha,ierr
real(rk)                       :: gm,a,e,inc,capom,omega,capm
integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after

c...  Executable code 

c...  if first time through open file
if (i1st.eq.0) then
   call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
   if (ierr.ne.0) then
      write(*,*) ' SWIFT ERROR: in io_write_frame: '
      write(*,*) '     Could not open binary output file:'
      call util_exit(1)
   endif
   i1st = 1
else
   call io_open(iu,oname,'append','UNFORMATTED',ierr)
endif

call io_write_hdr(iu,time,nbod,ntp,istat)
      
c...  write out planets
do i=2,nbod
   gm = mass(1)+mass(i)
   id = -1*i
   call orbel_xv2el(xh(1:3,i),vxh(1:3,i),gm,ialpha,a,e,inc,capom,omega,capm)
   call io_write_line(iu,id,a,e,inc,capom,omega,capm)
enddo

c...  write out test particles
gm = mass(1)
do i=1,ntp
   if (istat(i,1).eq.0) then
      call orbel_xv2el(xht(1:3,i),vxht(1:3,i),gm,ialpha,a,e,inc,capom,omega,capm)
      call io_write_line(iu,i,a,e,inc,capom,omega,capm)
   endif
enddo

close(iu)
return
end subroutine io_write_frame
