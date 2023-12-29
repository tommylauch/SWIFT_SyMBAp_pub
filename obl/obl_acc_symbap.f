c***************************************************************************
c                  OBL_ACC_SYMBAP.F
c*************************************************************************
c OBL_ACC returns the BARYCENTRIC x,y,z components of the acc. on NBOD
c particles due to the oblateness of mass(1) using  
c the values of J2RP2 and J4RP4 passed into the routine.
c (J2RP2 for example is the product of 
c J_2 times the square of the central body's radius)
c Here we return the net acc. produced
c only by the J2 and J4 terms (i.e. including
c neither the monopole nor higher order terms).
c      
c
c             Input:
c                 nbod     ==>  number of massive bodies (incl. central one)
c                 mass(*)  ==>  masses of particles (real*8 array)
c                 j2rp2    ==>  scaled value of j2 moment (real*8 scalar)
c                 j4rp4    ==>  scaled value of j4 moment (real*8 scalar)
c                                    (real*8 vectors)
c                 xh(*)    ==>  HELIO. positions of particles
c                 irh(*)   ==> 1./ magnitude of radius vector (real*8 vector)
c                                (passed in to save calcs.)
c             Output:
c               aoblx(*)   ==>  BARY. components of accel 
c                                        (real*8 vectors) 
c
c Remarks:  aoblx(1) (for example) contains x-component of
c           bary. acc. of central body
c Authors:  Martin Duncan 
c Date:    3/4/94
c Last revision: 

      subroutine obl_acc_symbap(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 j2rp2,j4rp4
      real*8 mass(NPLMAX)
      real*8 xh(3,NPLMAX),irh(NPLMAX)

c...  Output
      real*8 aoblx(3,NPLMAX)

c...  Internals
      integer n
      real*8 rinv2,t0,t1,t2,t3
      real*8 fac1,fac2

c----
c...  executable code

c First get the bary acc. of each "planet" due to central oblate "sun"
      do n=2,nbod
c Note that here we assume we know inverse of radius rather than calc. it
c from (x,y,z) to save the sqrt.
         rinv2 = irh(n)**2
         t0 = -mass(1)*rinv2*rinv2*irh(n)
         t1 = 1.5d0 *j2rp2
         t2 = xh(3,n)**2*rinv2
         t3 = 1.875d0 *j4rp4*rinv2

         fac1 = t0*(t1 - t3 - (5.d0*t1 - (14.d0 - 21.d0*t2)*t3)*t2)
         fac2 = 2.d0*t0*(t1 - (2.d0 - (14.d0*t2/3.d0))*t3)
      
         aoblx(1,n) = fac1*xh(1,n)
         aoblx(2,n) = fac1*xh(2,n)
         aoblx(3,n) = (fac1 + fac2)*xh(3,n)
      enddo
c Now compute the bary. acc. of Sun due to all the planets
      aoblx(1,1) = 0.d0
      aoblx(2,1) = 0.d0
      aoblx(3,1) = 0.d0
      do n=2,nbod
         aoblx(1,1) = aoblx(1,1) - mass(n)*aoblx(1,n)/mass(1)
         aoblx(2,1) = aoblx(2,1) - mass(n)*aoblx(2,n)/mass(1)
         aoblx(3,1) = aoblx(3,1) - mass(n)*aoblx(3,n)/mass(1)
      enddo

      return
      end                       !  obl_acc_symbap.f
c____________________________________________________________________________
