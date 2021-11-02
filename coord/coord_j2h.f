c***********************************************************************
c	                    COORD_J2H.F
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

	subroutine coord_j2h(nbod,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xh,yh,zh,vxh,vyh,vzh)


      include '../swift.inc'

c...  Inputs: 
	integer nbod 
	real*8 mass(nbod)
	real*8 xj(nbod),yj(nbod),zj(nbod)
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Outputs:
	real*8 xh(nbod),yh(nbod),zh(nbod)
	real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
	integer n
	real*8 sumx,sumy,sumz,sumvx,sumvy,sumvz
	real*8 eta(NTPMAX)

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi coords

	eta(1) = mass(1)
	do n = 2,nbod
	   eta(n) = eta(n-1) + mass(n)
        enddo

        xh(1) = 0.d0
        yh(1) =  0.d0
	zh(1) =  0.d0
	vxh(1) =  0.d0
	vyh(1) =  0.d0
	vzh(1) =  0.d0

        xh(2) = xj(2) 
        yh(2) = yj(2)
	zh(2) = zj(2)
	vxh(2) = vxj(2)
	vyh(2) = vyj(2)
	vzh(2) = vzj(2)

	sumx = mass(2)*xj(2)/eta(2)
	sumy = mass(2)*yj(2)/eta(2)
	sumz = mass(2)*zj(2)/eta(2)
	sumvx = mass(2)*vxj(2)/eta(2)
	sumvy = mass(2)*vyj(2)/eta(2)
	sumvz = mass(2)*vzj(2)/eta(2)

	do n=3,nbod 
	  xh(n) = xj(n) + sumx
	  yh(n) = yj(n) + sumy
	  zh(n) = zj(n) + sumz
	  vxh(n) = vxj(n) + sumvx
	  vyh(n) = vyj(n) + sumvy
	  vzh(n) = vzj(n) + sumvz

	  if(n.lt.nbod) then
             sumx = sumx + mass(n)*xj(n)/eta(n)
             sumy = sumy + mass(n)*yj(n)/eta(n)
             sumz = sumz + mass(n)*zj(n)/eta(n)
             sumvx = sumvx + mass(n)*vxj(n)/eta(n)
             sumvy = sumvy + mass(n)*vyj(n)/eta(n)
             sumvz = sumvz + mass(n)*vzj(n)/eta(n)
          endif

	enddo


123	return
	end     ! coord_j2h

c--------------------------------------------------------------------------

