!***************************************************************************
!                  OBL_ACC_SYMBAP.F
!*************************************************************************
! OBL_ACC returns the BARYCENTRIC x,y,z components of the acc. on NBOD
! particles due to the oblateness of mass(1) using  
! the values of J2RP2 and J4RP4 passed into the routine.
! (J2RP2 for example is the product of 
! J_2 times the square of the central body's radius)
! Here we return the net acc. produced
! only by the J2 and J4 terms (i.e. including
! neither the monopole nor higher order terms).
!             Input:
!                 nbod     ==>  number of massive bodies (incl. central one)
!                 mass(*)  ==>  masses of particles (real*8 array)
!                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
!                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
!                                    (real*8 vectors)
!                 xh(*)    ==>  HELIO. positions of particles
!                 irh(*)   ==> 1./ magnitude of radius vector (real*8 vector)
!                                (passed in to save calcs.)
!             Output:
!               aoblx(*)   ==>  BARY. components of accel 
!                                        (real*8 vectors) 
! Remarks:  aoblx(1) (for example) contains x-component of
!           bary. acc. of central body
! Authors:  Martin Duncan 
! Date:    3/4/94
! Last revision: 

subroutine obl_acc(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: j2rp2,j4rp4
real(rk), intent(in)    :: mass(:),xh(:,:),irh(:)

real(rk), intent(out)   :: aoblx(:,:)

integer(ik)             :: n
real(rk)                :: rinv2,t0,t1,t2,t3
real(rk)                :: fac1,fac2

!...  Executable code

! First get the bary acc. of each "planet" due to central oblate "sun"
   do n=2,nbod
! Note that here we assume we know inverse of radius rather than calc. it
! from (x,y,z) to save the sqrt.
      rinv2 = irh(n)**2
      t0 = -mass(1)*rinv2*rinv2*irh(n)
      t1 = 1.5_rk *j2rp2
      t2 = xh(3,n)**2*rinv2
      t3 = 1.875_rk *j4rp4*rinv2

      fac1 = t0*(t1-t3-(5.0_rk*t1 - (14.0_rk - 21.0_rk*t2)*t3)*t2)
      fac2 = 2.0_rk*t0*(t1 - (2.0_rk - (14.0_rk*t2/3.0_rk))*t3)

      aoblx(1,n) = fac1*xh(1,n)
      aoblx(2,n) = fac1*xh(2,n)
      aoblx(3,n) = (fac1 + fac2)*xh(3,n)
   enddo
! Now compute the bary. acc. of Sun due to all the planets
   aoblx(1,1) = 0.0_rk
   aoblx(2,1) = 0.0_rk
   aoblx(3,1) = 0.0_rk
   do n=2,nbod
      aoblx(:,1) = aoblx(:,1) - mass(n)*aoblx(:,n)/mass(1)
   enddo

return
end subroutine obl_acc
