!*************************************************************************
!                          ANAL_ENERGY_MTINY_SYMBAP.F
!*************************************************************************
! Calculates the energy of the total system (massive bodies) wrt time.
! returns the total energy of n objects by direct pairwise summation
! G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
!      Input:
!            t             ==>  current time
!            nbod          ==>  number of massive bodies (int scalar)
!            nbodm         ==>  Location of last massive body(int scalar)
!            mass          ==>  mass of bodies (real array)
!            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
!            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
!            xh            ==>  current position in heliocentric coord 
!                               (real array)
!            vxh           ==>  current velocity in heliocentric coord 
!                               (real array)
!      Output:
!            ke            ==>  kinetic energy
!            pot           ==>  potential energy
!            energy        ==>  Total energy
!            eltot         ==>  components of total angular momentum
!                               (real array)
! Remarks: Based on anal_energy
! Authors:  Hal Levison
! Date:  12/16/06
! Last revision:  

subroutine anal_energy_mtiny(nbod,nbodm,mass,j2rp2,j4rp4,xh,vxh,       &
                             ke,pot,energy,eltot)
use swift_mod
use coord_interface
use mvs_interface
implicit none

integer(ik), intent(in) :: nbod,nbodm
real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:),j2rp2,j4rp4

real(rk), intent(out)   :: energy,eltot(:),ke,pot

real(rk)                :: elx(3),xx(3)
real(rk)                :: rr2,oblpot,msys,irh(NTPMAX),ir3h(NTPMAX)
real(rk)                :: xb(3,NTPMAX),vxb(3,NTPMAX)
integer(ik)             :: i,j

!...  Executable code 

   call coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)   

   eltot(1) = (xb(2,nbod)*vxb(3,nbod)-xb(3,nbod)*vxb(2,nbod))
   eltot(2) = (xb(3,nbod)*vxb(1,nbod)-xb(1,nbod)*vxb(3,nbod))
   eltot(3) = (xb(1,nbod)*vxb(2,nbod)-xb(2,nbod)*vxb(1,nbod))
   eltot = eltot*mass(nbod)

   ke = 0.5_rk*mass(nbod)*dot_product(vxb(1:3,nbod),vxb(1:3,nbod))
   pot = 0.0_rk

   do i=1,nbodm
      elx(1) = xb(2,i)*vxb(3,i)-xb(3,i)*vxb(2,i)
      elx(2) = xb(3,i)*vxb(1,i)-xb(1,i)*vxb(3,i)
      elx(3) = xb(1,i)*vxb(2,i)-xb(2,i)*vxb(1,i)
      elx = elx*mass(i)
      eltot(:) = eltot(:) + elx(:)
      
      ke = ke+0.5_rk*mass(i)*dot_product(vxb(1:3,i),vxb(1:3,i))
      do j = i+1,nbod
         xx(:) = xb(:,i)-xb(:,j)
         rr2 = dot_product(xx,xx)
         if ((mass(i).ne.0.0_rk).and.(mass(j).ne.0.0_rk)) then
            pot = pot-mass(i)*mass(j)/(sqrt(rr2))
         endif
      enddo
   enddo

   do i=nbodm+1,nbod-1
      ke = ke + 0.5_rk*mass(i)*dot_product(vxb(1:3,i),vxb(1:3,i))
      elx(1) = xb(2,i)*vxb(3,i)-xb(3,i)*vxb(2,i)
      elx(2) = xb(3,i)*vxb(1,i)-xb(1,i)*vxb(3,i)
      elx(3) = xb(1,i)*vxb(2,i)-xb(2,i)*vxb(1,i)
      elx = elx*mass(i)
      eltot(:) = eltot(:) + elx(:)
   enddo

   if (j2rp2.ne.0.0_rk) then
      call getacch_ir3(nbod,2,xh,ir3h,irh)
      call obl_pot(nbod,mass,j2rp2,j4rp4,xh,irh,oblpot)
      pot = pot + oblpot
   endif

   energy = ke+pot

return
end subroutine anal_energy_mtiny
