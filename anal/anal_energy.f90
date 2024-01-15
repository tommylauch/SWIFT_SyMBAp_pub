c*************************************************************************
c                          ANAL_ENERGY_SYMBAP.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh            ==>  current position in heliocentric coord 
c                               (real array)
c            vxh           ==>  current velocity in heliocentric coord 
c                               (real array)
c
c      Output:
c            ke            ==>  kinetic energy
c            pot           ==>  potential energy
c            energy        ==>  Total energy
c            eltot         ==>  components of total angular momentum
c                               (real array)
c
c Remarks: 
c Authors:  Martin Duncan
c Date:  ?
c Last revision:  1/24/97 HFL 

subroutine anal_energy(nbod,mass,j2rp2,j4rp4,xh,vxh,ke,pot,energy,eltot)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:),j2rp2,j4rp4
real(rk), intent(in)    :: xh(:,:),vxh(:,:)

real(rk), intent(out)   :: energy,eltot(:),ke,pot

real(rk)                :: elx(3),xx(3)
real(rk)                :: xb(3,NTPMAX),vxb(3,NTPMAX)
real(rk)                :: rr2,oblpot,msys,irh(NTPMAX),ir3h(NTPMAX)
integer(ik)             :: i,j

c...  Executable code 

call coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)   

eltot(1) = xb(2,nbod)*vxb(3,nbod)-xb(3,nbod)*vxb(2,nbod)
eltot(2) = xb(3,nbod)*vxb(1,nbod)-xb(1,nbod)*vxb(3,nbod)
eltot(3) = xb(1,nbod)*vxb(2,nbod)-xb(2,nbod)*vxb(1,nbod)
eltot = eltot*mass(nbod)
      
ke = 0.5_rk*mass(nbod)*dot_product(vxb(1:3,nbod),vxb(1:3,nbod))
pot = 0.0_rk

do i=1,nbod-1
   elx(1) = xb(2,i)*vxb(3,i)-xb(3,i)*vxb(2,i)
   elx(2) = xb(3,i)*vxb(1,i)-xb(1,i)*vxb(3,i)
   elx(3) = xb(1,i)*vxb(2,i)-xb(2,i)*vxb(1,i)
   elx = elx*mass(i)
   eltot(:) = eltot(:) + elx(:)
   
   ke = ke+0.5_rk*mass(i)*dot_product(vxb(1:3,i),vxb(1:3,i))
   do j=i+1,nbod
      xx(:) = xb(:,i)-xb(:,j)
      rr2 = dot_product(xx,xx)
      if ( (mass(i).ne.0.0_rk).and.(mass(j).ne.0.0_rk) ) then
         pot = pot-mass(i)*mass(j)/(sqrt(rr2))
      endif
   enddo
enddo

if (j2rp2.ne.0.0_rk) then
   call getacch_ir3(nbod,2_ik,xh,ir3h,irh)
   call obl_pot(nbod,mass,j2rp2,j4rp4,xh,irh,oblpot)
   pot = pot+oblpot
endif

energy = ke+pot

return
end subroutine anal_energy
