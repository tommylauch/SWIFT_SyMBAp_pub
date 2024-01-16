!*************************************************************************
!                            DISCARD_MASS_MERGE5P_MTINY.F
!*************************************************************************
! Merge two massive bodies
!             Input:
!                 time          ==>  current time (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 ip1,ip2       ==>  planets to merge (real scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>   position in helio coord 
!                                    (real arrays)
!                 vxh           ==>   pl vel in helio coord 
!                                    (real arrays)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 eoff          ==> Amount of energy lost due to discards
!                                          (real scalar)
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
!             Output:
!                 mass          ==>  recalculated mass of bodies (real array)
!                 xh            ==>  recalculated position in helio coord 
!                                    (real arrays)
!                 vxh           ==>  recalculated pl vel in helio coord 
!                                    (real arrays)
!                 rpl           ==>  recalculated physical sizes of a planet.
!                                    (real array)
!                 eoff          ==> Updated amount of energy lost from discards
!                                          (real scalar)
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
! Remarks: Based on discard_mass_merge5
! Authors:  Hal Levison 
! Date:    12/16/09
! Last revision:  01/09/22 energy offset ignored for SyMBAp

subroutine discard_mass_merge5p_mtiny(time,nbod,nbodm,ip1,ip2,         &
           mass,xh,vxh,rpl,eoff,ielc,ielst)
implicit none
use swift_mod
use io_interface

integer(ik), intent(in)    :: ip1,ip2
real(rk), intent(in)       :: time

integer(ik), intent(inout) :: nbod,nbodm
real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:),rpl(:),eoff
integer(ik), intent(inout) :: ielst(:,:),ielc

real(rk)                   :: mtot,m(2),r(2),x(3,2),vx(3,2)
integer(ik)                :: itmp,j,i
real(rk)                   :: j2rp2,j4rp4

!...  Executable code

   j2rp2 = 0.0_rk
   j4rp4 = 0.0_rk

   if (mass(ip2).gt.mass(ip1)) then
      itmp = ip1
      ip1 = ip2
      ip2 = itmp
   endif

   write(*,*) ' Merging particles ',ip1, ' and ', ip2,' at t= ',time

   x(:,1) = xh(:,ip1)
   vx(:,1) = vxh(:,ip1)
   m(1) = mass(ip1)
   r(1) = rpl(ip1)

   x(:,2) = xh(:,ip2)
   vx(:,2) = vxh(:,ip2)
   m(2) = mass(ip2)
   r(2) = rpl(ip2)

!... Note:  I am just putting these guys together here, which is
!...        clearly wrong.  I should integrate back to the time
!...        of close approach.

   mtot = mass(ip1)+mass(ip2)

   rpl(ip1) = (r(1)**3 + r(2)**3)**ONETHRD
   vxh(:,ip1) = (mass(ip1)*vx(:,1)+mass(ip2)*vx(:,2))/mtot
   mass(ip1) = mtot

!..   Put in zeros for the rest the second particle
   xh(:,ip2) = xh(:,ip2)*1.0e10_rk   ! so danby does not fail
   vxh(:,ip2) = 0.0_rk
   mass(ip2) = 0.0_rk
   rpl(ip2) = 0.0_rk

!..   Remove any encounters with ip2
   j = 1_ik
   do while (j.le.ielc)
      if( (ielst(1,j).eq.ip2) .or. (ielst(2,j).eq.ip2) ) then
         do i=j+1,ielc
            ielst(1,i-1) = ielst(1,i)
            ielst(2,i-1) = ielst(2,i)
         enddo
         ielc = ielc - 1
      else
         j = j + 1
      endif
   enddo

   eoff = 0.0_rk

   call io_discard_merge(time,ip1,ip2,m,r,x,vx,                        &
        mass(ip1),rpl(ip1),xh(1:3,ip1),vxh(1:3,ip1))

return
end subroutine discard_mass_merge5p_mtiny
