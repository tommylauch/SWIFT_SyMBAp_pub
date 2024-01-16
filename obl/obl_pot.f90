!***************************************************************************
!                  OBL_POT_SYMBAP.F
!*************************************************************************
! OBL_POT returns the total potential in the barycentric frame for NBOD
! particles due to the oblateness of mass(1) using  
! the values of J2RP2 and J4RP4 passed into the routine.
! (J2RP2 for example is the product of 
! J_2 times the square of the central body's radius)
! Here we return the potential produced
! only by the J2 and J4 terms (i.e. including
! neither the monopole nor higher order terms).
!             Input:
!                 nbod     ==>  number of massive bodies (incl. central one)
!                 mass(*)  ==>  masses of particles (real*8 array)
!                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
!                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
!                 xh(3,*)   ==>  HELIO. positions of particles
!                                    (real*8 vector)
!                 irh(*)   ==> 1./ magnitude of radius vector (real*8 vector)
!                                (passed in to save calcs.)
!             Output:
!                 oblpot  ==>  BARY. potential
!                                        (real*8 scalar) 
! Remarks:  
! Authors:  Martin Duncan 
! Date:    3/4/94
! Last revision: 

subroutine obl_pot(nbod,mass,j2rp2,j4rp4,xh,irh,oblpot)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: mass(:)
real(rk), intent(in)    :: j2rp2,j4rp4
real(rk), intent(in)    :: xh(:,:),irh(:)

real(rk), intent(out)   :: oblpot

integer(ik)             :: n
real(rk)                :: rinv2,t0,t1,t2,t3,p2,p4

!...  Executable code
! Sum all the the bary terms for each "planet" due to central oblate "sun"
   oblpot = 0.0_rk
   do n=2,nbod
! Note that here we assume we know inverse of radius rather than calc. it
! from (x,y,z) to save the sqrt.
      rinv2 = irh(n)**2
      t0 = mass(1)*mass(n)*rinv2*irh(n)
      t1 = j2rp2
      t2 = xh(3,n)**2*rinv2
      t3 = j4rp4*rinv2

      p2 = 0.5_rk*(3.0_rk*t2-1.0_rk)
      p4 = 0.125_rk*((35.0_rk*t2-30.0_rk)*t2+3.0_rk)

      oblpot = oblpot+t0*(t1*p2+t3*p4)
   enddo

return
end subroutine obl_pot

