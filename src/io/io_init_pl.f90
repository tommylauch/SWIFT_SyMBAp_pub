!************************************************************************
!                             IO_INIT_PL.F
!************************************************************************
! IO_INIT_PL reads in the data for the Sun and planets for 
! symba routines
!            Input:
!                infile        ==> File name to read from (character*80)
!                lclose        ==> .true. --> discard particle if it gets 
!                                   too close to a planet. Read in that 
!                                   distance in io_init_pl_symba
!                                     (logical*2 scalar)
!                iflgchk        ==>  bit 5 set ==>  include J2 and J4 terms
!            Output:
!                nbod          ==>  number of massive bodies (int scalar)
!                mass          ==>  mass of bodies (real array)
!                xh            ==>  initial position in Helio coord 
!                                   (real array)
!                vxh           ==>  initial position in Helio coord 
!                                   (real array)
!                rpl           ==>  physical size of planet
!                                   (real array)
!                rhill         ==>  size of planet's hills sphere
!                                   (real array)
!                j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                    (real scalars)
! Remarks: Based on io_init_pl 
! Authors:  Hal Levison
! Date:    11/21/96
! Last revision: 1/10/97

subroutine io_init_pl(infile,lclose,iflgchk,nbod,mass,xh,vxh,rpl,rhill,&
                      j2rp2,j4rp4)
use swift_mod
use io_mod
use util_interface
use io_interface, except_this_one => io_init_pl
implicit none

character(len=*), intent(in)   :: infile
integer(ik), intent(in)        :: iflgchk
logical(ik), intent(in)        :: lclose

real(rk), intent(out)          :: mass(:),rpl(:),j2rp2,j4rp4
real(rk), intent(out)          :: xh(:,:),vxh(:,:),rhill(:)
integer(ik), intent(out)       :: nbod

integer(ik)                    :: j,ierr,ibad
real(rk)                       :: r2hill(NTPMAX),rhrat

!...  Executable code

   write(*,*) 'Planet data file is ',infile
   call io_open(7,infile,'old','formatted',ierr)

! Read number of planets
   read(7,*) nbod

   if (nbod .gt. NTPMAX) then
      write(*,*) ' SWIFT ERROR: in io_init_pl_symbap: '
      write(*,*) '   The number of massive bodies,',nbod,','
      write(*,*) '   is too large, it must be less than',NTPMAX
      call util_exit(1)
   endif

   write(*,*) 'Number of bodies (incl. the Sun) is ',nbod

! For each planet read mass, and helioc. position and vel .
   if (btest(iflgchk,5)) then ! bit 5 is set
      read(7,*) mass(1),j2rp2,j4rp4
   else
      read(7,*) mass(1)
      j2rp2 = 0.0d0
      j4rp4 = 0.0d0
   endif
   read(7,*) xh(1,1),xh(2,1),xh(3,1)
   read(7,*) vxh(1,1),vxh(2,1),vxh(3,1)
   rpl(1) = 0.0d0
   rhill(1) = 0.0d0

   if ( (xh(1,1).ne.0.0d0) .or.                                        &
       (xh(2,1).ne.0.0d0) .or.                                         &
       (xh(3,1).ne.0.0d0) .or.                                         &
       (vxh(1,1).ne.0.0d0) .or.                                        &
       (vxh(2,1).ne.0.0d0) .or.                                        &
       (vxh(3,1).ne.0.0d0) ) then
      write(*,*) 'SWIFT ERROR: in io_init_pl_symbap: '
      write(*,*) '  Input MUST be in heliocentric coordinates '
      write(*,*) '  Position and Vel. of Massive body 1 .ne. 0'
      call util_exit(1)
   endif

   do j=2,nbod
      if (lclose) then
         read(7,*) mass(j),rhill(j),rpl(j)
      else
         read(7,*) mass(j),rhill(j)
      endif
      read(7,*) xh(1,j),xh(2,j),xh(3,j)
      read(7,*) vxh(1,j),vxh(2,j),vxh(3,j)
   enddo

close(unit = 7)

!...  check to see if the hills spheres are ok
   call util_hills(nbod,mass,xh,vxh,r2hill) 
   ibad = 0
   do j=2,nbod
      rhrat = rhill(j)/sqrt(r2hill(j))
      if ( (rhrat.gt.2.0) .or. (rhrat.lt.0.5) ) ibad = ibad + 1
   enddo

   if (ibad.ne.0) then
      write(*,*) 'Warning in io_init_pl_symbap:'
      write(*,*) '  Hill''s spheres are not consistent on ',ibad,' objects'
   endif
      
return
end subroutine io_init_pl
