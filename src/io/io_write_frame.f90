!*************************************************************************
!                           IO_WRITE_FRAME
!*************************************************************************
! write out a whole frame to an real*4 binary file.
! both massive and test particles
!            Input:
!                time          ==>  current time (real scalar)
!                nbod          ==>  number of massive bodies (int scalar)
!                ntp           ==>  number of massive bodies (int scalar)
!                mass          ==>  mass of bodies (real array)
!                xh            ==>  current position in helio coord 
!                                   (real array)
!                vxh           ==>  current velocity in helio coord 
!                                   (real array)
!                xht           ==>  current part position in helio coord 
!                                   (real array)
!                vxht          ==>  current velocity in helio coord 
!                                   (real array)
!                istat         ==>  status of the test paricles
!                                   (2d integer array)
!                                   istat(i,1) = 0 ==> active:  = 1 not
!                                   istat(i,2) = -1 ==> Danby did not work
!                oname         ==>  output file name (character string) 
!                iu            ==>  unit number to write to
!                fopenstat     ==>  The status flag for the open 
!                                   statements of the output files.  
!                                   (character*80)
! Remarks: Based on io_write_frame
! Authors:  Hal Levison 
! Date:    2/22/94
! Last revision: 

subroutine io_write_frame(time,nbod,ntp,mass,xh,vxh,xht,vxht,istat,    &
                          oname,iu,fopenstat)
use swift_mod
use util_interface
use io_interface, except_this_one => io_write_frame
use orbel_interface
implicit none

integer(ik), intent(in)        :: nbod,ntp,iu
real(rk), intent(in)           :: mass(:),time
integer(ik), intent(in)        :: istat(:,:)
real(rk), intent(in)           :: xh(:,:),vxh(:,:)
real(rk), intent(in)           :: xht(:,:),vxht(:,:)
character(len=50), intent(in)  :: oname,fopenstat

integer(ik)                    :: i,id
integer(ik)                    :: ialpha,ierr
real(rk)                       :: gm,a,e,inc,capom,omega,capm
integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after

!...  Executable code

!...  if first time through open file
   if (i1st.eq.0) then
      call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
      if (ierr.ne.0) then
         write(*,*) ' SWIFT ERROR: in io_write_frame: '
         write(*,*) '     Could not open binary output file:'
         call util_exit(1)
      endif
      i1st = 1_ik
   else
      call io_open(iu,oname,'append','UNFORMATTED',ierr)
   endif

   call io_write_hdr(iu,time,nbod,ntp,istat)
      
!...  write out planets
   do i=2,nbod
      gm = mass(1)+mass(i)
      id = -1*i
      call orbel_xv2el(xh(1:3,i),vxh(1:3,i),gm,ialpha,a,e,inc,capom,omega,capm)
      call io_write_line(iu,id,a,e,inc,capom,omega,capm)
   enddo

!...  write out test particles
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
