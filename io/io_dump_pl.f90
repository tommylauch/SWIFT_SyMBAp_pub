!***********************************************************************
!                         IO_DUMP_PL.F
!***********************************************************************
! Dumps the data for the Sun and planets 
!             Input:
!                 dplfile       ==>  Name of file to write to (character*80)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>  initial position in Helio coord 
!                                    (real arrays)
!                 vxh           ==>  initial position in Helio coord 
!                                    (real arrays)
!                 lclose        ==> .true. --> discard particle if it gets 
!                                    too close to a planet. Read in that 
!                                    distance in io_init_pl
!                                      (logical*2 scalar)
!                 iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
!                 rpl           ==>  physical size of planet
!                                    (real array)
!                 rhill         ==>  size of planet's hills sphere
!                                    (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
! Remarks: Based on io_dump_pl.f
! Authors:  Hal Levison
! Date:    1/8/97
! Last revision: 

subroutine io_dump_pl(dplfile,nbod,mass,xh,vxh,lclose,                 &
                      iflgchk,rpl,rhill,j2rp2,j4rp4)
implicit none
use swift_mod
use io_interface, except_this_one => io_dump_pl

integer(ik), intent(in)        :: nbod,iflgchk
real(rk), intent(in)           :: mass(:),rpl(:),j2rp2,j4rp4
real(rk), intent(in)           :: xh(:,:),vxh(:,:),rhill(:)
character(len = :), intent(in) :: dplfile
logical(ik), intent(in)        :: lclose

integer(ik)                    :: j,ierr

!...  Executable code

   call io_open(7,dplfile,'unknown','formatted',ierr)

   write(7,*) nbod

   if (btest(iflgchk,5)) then ! bit 5 is set
      write(7,fmt_rrr) mass(1),j2rp2,j4rp4
   else
      write(7,fmt_rrr) mass(1)
   endif

   write(7,fmt_rrr) xh(1,1),xh(2,1),xh(3,1)
   write(7,fmt_rrr) vxh(1,1),vxh(2,1),vxh(3,1)

   do j=2,nbod
      if (lclose) then
         write(7,fmt_rrr) mass(j),rhill(j),rpl(j)
      else
         write(7,fmt_rrr) mass(j),rhill(j)
      endif
      write(7,fmt_rrr) xh(1,j),xh(2,j),xh(3,j)
      write(7,fmt_rrr) vxh(1,j),vxh(2,j),vxh(3,j)
   enddo

   close(unit = 7)

return
end subroutine io_dump_pl
