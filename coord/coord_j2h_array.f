c***********************************************************************
c	                    COORD_J2H_ARRAY.F
c***********************************************************************
*     PURPOSE: Converts from Jacobi to Helio coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xj(*),yj(*),zj(*) ==> Jacobi particle coords
*                                          (real array)
*		     vxj(*),vyj(*),vzj(*) ==> Jacobi particle velocities
*                                             (real array)
*                 Returned are
*                    xh(*),yh(*),zh(*) ==> Helio particle positions
*                                          (real array)
*                    vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                            (real array)
*       
*     ALGORITHM: See my notes on Nov 21. 
*
*     Authors:  Martin Duncan
*     WRITTEN:  Jan 27/93
*     REVISIONS: 2/20/2K  HFL

	subroutine coord_j2h_array(nbod,mass,xj,vxj,xh,vxh)


      include '../swift.inc'

c...  Inputs: 
      integer nbod 
      real*8 mass(nbod),xj(3,nbod),vxj(3,nbod)

c...  Outputs:
      real*8 xh(3,nbod),vxh(3,nbod)

c...  Internals:
      integer n
      real*8 sumx(3),sumvx(3)
      real*8 eta(NTPMAX)

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

      eta(1) = mass(1)
      do n = 2,nbod
         eta(n) = eta(n-1) + mass(n)
      enddo

      xh(:,1) = 0.d0
      vxh(:,1) =  0.d0

      xh(:,2) = xj(:,2) 
      vxh(:,2) = vxj(:,2)

      sumx = mass(2)*xj(:,2)/eta(2)
      sumvx = mass(2)*vxj(:,2)/eta(2)

      do n=3,nbod 
         xh(:,n) = xj(:,n) + sumx
         vxh(:,n) = vxj(:,n) + sumvx

         if (n.lt.nbod) then
            sumx = sumx + mass(n)*xj(:,n)/eta(n)
            sumvx = sumvx + mass(n)*vxj(:,n)/eta(n)
         endif
      enddo
      return
      end     ! coord_j2h_array

c--------------------------------------------------------------------------

