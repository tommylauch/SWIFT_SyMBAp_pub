***********************************************************************
c                      ORBEL_XV2AEQ.F
***********************************************************************
*     PURPOSE:  Given the cartesian position and velocity of an orbit,
*       compute the osculating orbital elements a, e, and q only.
*
C       input:
c            x        ==>  position of object (real scalars)
c            vx       ==>  velocity of object (real scalars)
c            gmsum    ==> G*(M1+M2) (real scalar)
c
c       Output:
c           ialpha    ==> conic section type ( see PURPOSE, integer scalar)
C           a         ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
c            q        ==> perihelion distance (real scalar); q = a(1 - e)
c
*     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
*     REMARKS: Based on M. Duncan's orbel_xv2el.f
*      This routine is generally applied to study (hyperbolic) close 
c       encounters of test particles with planets.
*     AUTHOR:  L. Dones
*     DATE WRITTEN:  February 24, 1994
*     REVISIONS: 
***********************************************************************

      subroutine orbel_xv2aeq_symbap(x,vx,gmsum,ialpha,a,e,q)

      include '../swift.inc'

c...  Inputs Only: 
      real*8 x(3),vx(3),gmsum

c...  Outputs
      integer ialpha
      real*8 a,e,q

c...  Internals:
      real*8 hx(3),h2,r,v2,energy,fac
c----
c...  Executable code 

* Compute the angular momentum H, and thereby the inclination INC.

      hx(1) = x(2)*vx(3) - x(3)*vx(2)
      hx(2) = x(3)*vx(1) - x(1)*vx(3)
      hx(3) = x(1)*vx(2) - x(2)*vx(1)
      h2 = hx(1)**2 + hx(2)**2 + hx(3)**2

*  Compute the radius R and velocity squared V2, and the dot
*  product RDOTV, the energy per unit mass ENERGY .

      r = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      v2 = vx(1)**2 + vx(2)**2 + vx(3)**2
      energy = 0.5d0*v2 - gmsum/r

*  Determine type of conic section and label it via IALPHA
      if(abs(energy*r/gmsum) .lt. sqrt(TINY)) then
         ialpha = 0
      else
         if(energy .lt. 0.d0) ialpha = -1 
         if(energy .gt. 0.d0) ialpha = +1
      endif

* Depending on the conic type, determine the remaining elements
***
c ELLIPSE :
      if(ialpha .eq. -1) then
        a = -0.5d0*gmsum/energy  
        fac = 1.d0 - h2/(gmsum*a)

         if (fac .gt. TINY) then
            e = sqrt ( fac )
        else
          e = 0.d0
        endif

          q = a*(1.d0 - e)
          endif
***
***
c HYPERBOLA
      if(ialpha .eq. +1) then
        a = +0.5d0*gmsum/energy  
        fac = h2/(gmsum*a)
          if (fac .gt. TINY) then
           e = sqrt ( 1.d0 + fac )
            q = -a*(1.d0 - e)
c     have to insert minus sign in expression for q because this code
c      takes a > 0, even for a hyperbola
        else
c we only get here if a hyperbola is essentially a parabola
c so we calculate e accordingly to avoid singularities
          e = 1.d0
            q = 0.5*h2/gmsum
        endif
          endif
***
***
c PARABOLA : ( NOTE - in this case "a", which is formally infinite,
c         is arbitrarily set equal to the pericentric distance q).
      if(ialpha .eq. 0) then
        a =  0.5d0*h2/gmsum  
        e = 1.d0
          q = a
      endif
***
***
      return
      end    ! orbel_xv2aeq
c------------------------------------------------------------------
