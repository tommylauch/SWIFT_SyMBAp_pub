***********************************************************************
*	                    COORD_H2J_ARRAY.F
***********************************************************************
*     PURPOSE: Converts from Heliocentric to Jacobi coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of  objects (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==> planetary masses (real array)
*		     xh(*),yh(*),zh(*) ==> barycentric particle coords
*                                            (real array)
*		     vxh(*),vyh(*),vzh(*) ==>barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    xj(*),yj(*),zj(*) ==> jacobi. particle positions
*                                             (real array)
*                    vxj(*),vyj(*),vzj(*) ==> jacobi. particle velocities
*                                              (real array)
*       
*     ALGORITHM: See my notes Nov 21/92 
*     REMARKS:  Note that we set the Jacobi coord of the Sun = 0
*               This is not in accord with the definition, but
*               since we never use the Sun's Jacobi coord, it was the
*               fastest thing to do.
*       
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  Jan 27, 1993.
*     REVISIONS:  2/20/2K HFL

	subroutine coord_h2j_array(nbod,mass,xh,vxh,xj,vxj)

      include '../swift.inc'

c...  Inputs: 
	integer nbod
	real*8 mass(nbod)
	real*8 xh(3,nbod),vxh(3,nbod)

c...  Outputs:
	real*8 xj(3,nbod),vxj(3,nbod)

c...  Internals:
	real*8 eta(NTPMAX)
	real*8 sumx(3),sumvx(3)
	real*8 capx(3),capvx(3)
	integer n

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

      eta(1) = mass(1)
      do n = 2,nbod
         eta(n) = eta(n-1) + mass(n)
      enddo

      xj(:,1) = 0.d0
      vxj(:,1) = 0.d0

      xj(:,2) = xh(:,2)
      vxj(:,2) = vxh(:,2)

      sumx = mass(2)*xh(:,2)
      sumvx = mass(2)*vxh(:,2)

      capx = sumx/eta(2)
      capvx = sumvx/eta(2)

      do n=3,nbod
         xj(:,n) = xh(:,n) - capx
         vxj(:,n) = vxh(:,n) - capvx

         if(n.lt.nbod) then
            sumx = sumx + mass(n)*xh(:,n)
            sumvx = sumvx + mass(n)*vxh(:,n)
            capx = sumx/eta(n)
            capvx = sumvx/eta(n)
         endif

      enddo

      return
	   end    ! coord_h2j_array
c--------------------------------------------------------------------
