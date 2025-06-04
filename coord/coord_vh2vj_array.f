***********************************************************************
*	                    COORD_VH2VJ_ARRAY.F
***********************************************************************
*     PURPOSE: Converts from Heliocentric to Jacobi coords. VELOCITIES ONLY.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of  objects (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==> planetary masses (real array)
*		     vxh(*),vyh(*),vzh(*) ==>barycentric particle velocities
*                                             (real array)
*                 Returned are
*                    vxj(*),vyj(*),vzj(*) ==> jacobi. particle velocities
*                                              (real array)
*       
*     ALGORITHM: See my notes Nov 21/92 
*     REMARKS:  
*       
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  Jan 29, 1993.
*     REVISIONS:  2/20/2K HFL

      subroutine coord_vh2vj_array(nbod,mass,vxh,vxj)


      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod)
      real*8 vxh(3,nbod)

c...  Outputs:
      real*8 vxj(3,nbod)

c...  Internals:
      real*8 eta(NTPMAX)
      real*8 sumvx(3)
      real*8 capvx(3)
      integer n

c----
c...  Executable code 

c First calc. the array eta(*) then convert to jacobi velocities

      eta(1) = mass(1)
      do n = 2,nbod
         eta(n) = eta(n-1) + mass(n)
      enddo

      vxj(:,1) = 0.d0
      vxj(:,2) = vxh(:,2)

      sumvx = mass(2)*vxh(:,2)
      capvx = sumvx/eta(2)

      do n=3,nbod
         vxj(:,n) = vxh(:,n) - capvx
         if(n.lt.nbod) then
            sumvx = sumvx + mass(n)*vxh(:,n)
            capvx = sumvx/eta(n)
         endif
      enddo

      return
      end    ! coord_vh2vj
c--------------------------------------------------------------------------
