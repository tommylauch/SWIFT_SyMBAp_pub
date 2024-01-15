c***********************************************************************
c 	               ORBEL_XV2AQT.F
c***********************************************************************
c     PURPOSE:  Given the cartesian position and velocity of an orbit,
c       compute the osculating orbital elements a, e, and q only.
c
c       input:
c            x,y,z    ==>  position of object (real scalars)
c            vx,vy,vz ==>  velocity of object (real scalars)
c            gmsum       ==> G*(M1+M2) (real scalar)
c
c       Output:
c	     ialpha   ==> conic section type ( see PURPOSE, integer scalar)
c	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            q        ==> perihelion distance (real scalar); q = a(1 - e)
c          capm       ==> mean anomoly(real scalar)
c          tperi      ==> time to next or last perihelion, which ever is less 
c                         (real scalar)
c
c     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
c     REMARKS: Based on M. Duncan's orbel_xv2el.f
c      This routine is generally applied to study (hyperbolic) close 
c       encounters of test particles with planets.
c     AUTHOR:  Hal Levison
c     DATE WRITTEN:  8/7/01
c     REMARKS: The tperi may not be correct for parabolic orbits. 
c              I Think it is OK but beware!
c     REVISIONS: 
c***********************************************************************

      subroutine orbel_xv2aqt(x,y,z,vx,vy,vz,gmsum,
     &     ialpha,a,q,capm,tperi)
      
      include '../swift.inc'
      
c...  Inputs Only: 
      real*8 x,y,z,vx,vy,vz,gmsum
      
c...  Outputs
      integer ialpha
      real*8 a,q,capm,tperi
      
c...  Internals:
      real*8 hx,hy,hz,h2,r,v2,energy,fac,vdotr,cape,e
      real*8 capf,tmpf,meanmo,face,w

c---- 
c...  Executable code 
      
c     Compute the angular momentum H, and thereby the inclination INC.
      
      hx = y*vz - z*vy
      hy = z*vx - x*vz
      hz = x*vy - y*vx
      h2 = hx*hx + hy*hy + hz*hz
      
*     Compute the radius R and velocity squared V2, and the dot
*     product RDOTV, the energy per unit mass ENERGY .
      
      r = sqrt(x*x + y*y + z*z)
      v2 = vx*vx + vy*vy + vz*vz
      energy = 0.5d0*v2 - gmsum/r
      vdotr = x*vx + y*vy + z*vz
      
*     Determine type of conic section and label it via IALPHA
      if(abs(energy*r/gmsum) .lt. sqrt(TINY)) then
         ialpha = 0
      else
         if(energy .lt. 0.d0) ialpha = -1 
         if(energy .gt. 0.d0) ialpha = +1
      endif
      
*     Depending on the conic type, determine the remaining elements
***   
c     ELLIPSE :
      if(ialpha .eq. -1) then
         a = -0.5d0*gmsum/energy  
         fac = 1.d0 - h2/(gmsum*a)
         
         if (fac .gt. TINY) then
            e = sqrt ( fac )
            face =(a-r)/(a*e)
            
c...  Apr. 16/93 : watch for case where face is slightly outside unity
            if ( face .gt. 1.d0) then
               cape = 0.d0
            else
               if ( face .gt. -1.d0) then
                  cape = acos( face )
               else
                  cape = PI
               endif
            endif
            If ( vdotr .lt. 0.d0 ) cape = 2.d0*PI - cape
         else
	    e = 0.d0
	    cape = 0.0d0
         endif
         
         capm = cape - e*sin (cape)
         q = a*(1.d0 - e)
      endif
***   
***   
c     HYPERBOLA
      if(ialpha .eq. +1) then
         a = +0.5d0*gmsum/energy  
         fac = h2/(gmsum*a)
         
         if (fac .gt. TINY) then
 	    e = sqrt ( 1.d0 + fac )
            q = -a*(1.d0 - e)
c     have to insert minus sign in expression for q because this code
c     takes a > 0, even for a hyperbola
            
	    tmpf = (a+r)/(a*e)
            if(tmpf.lt.1.0d0) then
               tmpf = 1.0d0
            endif
	    capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
	    if ( vdotr .lt. 0.d0 ) capf = - capf
            
         else
c     we only get here if a hyperbola is essentially a parabola
c     so we calculate e accordingly to avoid singularities
	    e = 1.d0
            q = 0.5*h2/gmsum
            
	    tmpf = (a+r)/(a*e)
	    capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
            
         endif
         
         capm = e * sinh(capf) - capf
                  
      endif
***
***
c PARABOLA : ( NOTE - in this case "a", which is formally infinite,
c         is arbitrarily set equal to the pericentric distance q).
      if(ialpha .eq. 0) then
         a =  0.5d0*h2/gmsum  
         e = 1.d0
         q = a
         w = acos(2.d0*a/r -1.d0)
         if ( vdotr .lt. 0.d0) w = 2.d0*PI - w
         tmpf = tan(0.5d0 * w)
         capm = tmpf* (1.d0 + tmpf*tmpf/3.d0)
      endif
***   
***   

      meanmo = sqrt(gmsum/(a**3))

      if( (capm.lt.PI) .or. (ialpha.ge.0) ) then
         tperi = -1.0d0*capm/meanmo
      else
         tperi = -1.0d0*(capm-TWOPI)/meanmo
      endif

      return
      end                       ! orbel_xv2aqt
c------------------------------------------------------------------
      
