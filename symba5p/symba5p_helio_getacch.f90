!*************************************************************************
!                        SYMBA5P_HELIO_GETACCH.F
!*************************************************************************
! This subroutine calculates the acceleration on the massive particles
! in the HELIOCENTRIC frame. 
!             Input:
!                 iflg        ==>  =0 calculate forces (int scalor)
!                                  =1 don't
!                 nbod        ==>  number of massive bodies (int scalor)
!                 nbodm       ==>  The last massive particle
!                                  (int scalor)
!                 mass        ==>  mass of bodies (real array)
!                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh          ==>  position in heliocentric coord (real arrays)
!             Output:
!                 axh         ==>  acceleration in helio coord (real arrays)
! Remarks Based on helio_getacch.f
! Author:  Hal Levison  
! Date:    9/12/99
! Last revision: 11/08/13 

subroutine symba5p_helio_getacch(iflg,nbod,nbodm,mass,j2rp2,j4rp4,     &
                                 xh,axh)
implicit none
use swift_mod
use sybam5p_mod
use mvs_interface
use obl_interface

integer(ik), intent(in) :: nbod,nbodm,iflg
real(rk), intent(in)    :: mass(:),xh(:,:),j2rp2,j4rp4

real(rk), intent(out)   :: axh(:,:)

integer(ik)             :: i,j
real(rk)                :: aoblx(3,NTPMAX)
real(rk), save          :: axhl(3,NTPMAX)
real(rk)                :: ir3h(NTPMAX),irh(NTPMAX)
real(rk)                :: dx(3),rji2,faci,facj,irij3

!...  Executable code

   if (iflg.eq.0) then
      axhl = 0.0_rk
!...     now the third terms
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& REDUCTION(+:axhl)
!$OMP& PRIVATE(i,j,dx,rji2,irij3,faci,facj)
!$OMP& SHARED(nbod,nbodm,mass,xh)
!$OMP DO COLLAPSE(2)
      do i=2,nbodm
         do j=i+1,nbod
            dx(:) = xh(:,j) - xh(:,i)
            rji2 = dot_product(dx,dx)

            irij3 = 1.0_rk/(rji2*sqrt(rji2))
            faci = mass(i)*irij3
            facj = mass(j)*irij3

            axhl(:,j) = axhl(:,j) - faci*dx(:)
            axhl(:,i) = axhl(:,i) + facj*dx(:)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif
!...  Now do j2 and j4 stuff
   if (j2rp2.ne.0.0_rk) then
      call getacch_ir3(nbod,2,xh,ir3h,irh)
      call obl_acc(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
      do i=2,nbod
         axh(:,i) = axhl(:,i) + aoblx(:,i)
      enddo
   else
      do i=2,nbod
         axh(:,i) = axhl(:,i)
      enddo
   endif

return
end subroutine symba5p_helio_getacch
