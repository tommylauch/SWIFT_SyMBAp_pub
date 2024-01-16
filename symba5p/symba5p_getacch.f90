!*************************************************************************
!                        SYMBA5P_GETACCH.F
!*************************************************************************
! This subroutine calculates the acceleration on the massive particles
! in the HELIOCENTRIC frame. 
!             Input:
!                 nbod        ==>  number of massive bodies (int scalor)
!                 nbodm       ==>  Location of last massive body(int scalar)
!                 mass        ==>  mass of bodies (real array)
!                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh          ==>  position in heliocentric coord (real arrays)
!                 mtiny       ==>  Small mass  (real array)
!                ielc         ==>  number of encounters (integer scalar)
!                ielst        ==>  list of ecnounters (2D integer array)
!             Output
!                 axh         ==>  acceleration in helio coord (real arrays)
! Remarks: Based on helio_getacch.f, but does not include the forces of
!          an bodxy B on body A, if body B and A are having an encounter.
! Author:  Hal Levison  
! Date:    3/20/97
! Last revision: 17/8/20

subroutine symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,axh,         &
                           mtiny,ielc,ielst)
implicit none
use swift_mod
use sybam5p_mod
use mvs_interface
use obl_interface

integer(ik), intent(in) :: nbod,nbodm,ielst(:,:),ielc
real(rk), intent(in)    :: mass(:),xh(:,:),j2rp2,j4rp4,mtiny

real(rk), intent(out)   :: axh(:,:)

real(rk)                :: aoblx(3,NTPMAX)
real(rk)                :: ir3h(NTPMAX),irh(NTPMAX) 
integer(ik)             :: i,j,ie
real(rk)                :: dx(3),rji2,faci,facj,irij3

!...  Executable code 

   axh = 0.0_rk
!...  now the third terms
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& REDUCTION(+:axh)
!$OMP& PRIVATE(i,j,ie,dx,rji2,irij3,faci,facj)
!$OMP& SHARED(nbod,nbodm,mass,xh,ielc,ielst)
!$OMP DO COLLAPSE(2)
   do i=2,nbodm
      do j=i+1,nbod
         dx(:) = xh(:,j) - xh(:,i)
         rji2 = dot_product(dx,dx)

         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3

         axh(:,j) = axh(:,j) - faci*dx(:)
         axh(:,i) = axh(:,i) + facj*dx(:)
      enddo
   enddo
!$OMP END DO NOWAIT
!...  Now subtract off anyone in an encounter
!$OMP DO
   do ie=1,ielc
      i = ielst(1,ie)
      j = ielst(2,ie)

      dx(:) = xh(:,j) - xh(:,i)
      rji2 = dot_product(dx,dx)

      irij3 = 1.0d0/(rji2*sqrt(rji2))
      faci = mass(i)*irij3
      facj = mass(j)*irij3

      axh(:,j) = axh(:,j) + faci*dx(:)
      axh(:,i) = axh(:,i) - facj*dx(:)
   enddo
!$OMP END DO
!$OMP END PARALLEL

!...  Now do j2 and j4 stuff
   if(j2rp2.ne.0.0_rk) then
      call getacch_ir3(nbod,2_ik,xh,ir3h,irh)
      call obl_acc(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
      do i=2,nbod
         if(mass(i).ne.0.0_rk) axh(:,i) = axh(:,i)+aoblx(:,i)
      enddo
   endif

return
end subroutine symba5p_getacch
