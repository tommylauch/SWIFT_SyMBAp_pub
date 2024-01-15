***********************************************************************
c                      ORBEL_XV2EL_SYMBAP.F
***********************************************************************
*     PURPOSE:  Given the cartesian position and velocity of an orbit,
*       compute the osculating orbital elements.
*
C       input:
c            x,y,z    ==>  position of object (real scalars)
c            vx,vy,vz ==>  velocity of object (real scalars)
c            gmsum       ==> G*(M1+M2) (real scalar)
c
c       Output:
c           ialpha   ==> conic section type ( see PURPOSE, integer scalar)
C           a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C           omega    ==> argument of perihelion (real scalar)
C           capm     ==> mean anomoly(real scalar)
c
*     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
*     REMARKS:  If the inclination INC is less than TINY, we
*       arbitrarily choose the longitude of the ascending node LGNODE
*       to be 0.0 (so the ascending node is then along the X axis).  If 
*       the  eccentricity E is less than SQRT(TINY), we arbitrarily
*       choose the argument of perihelion to be 0.
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 8,1992.
*     REVISIONS: 12/8/2011
***********************************************************************

      subroutine orbel_xv2el_symbap(x,vx,gmsum,ialpha,a,e,inc,capom,
     &                              omega,capm)
      include '../swift.inc'

c...  Inputs Only: 
      real*8 x(3),vx(3),gmsum

c...  Outputs
      integer ialpha
      real*8 a,e,inc,capom,omega,capm

c...  Internals:
      real*8 hx(3),h2,h,r,v2,v,vdotr,energy,fac,face,cape,capf,tmpf
      real*8 cw,sw,w,u

c----
c...  Executable code 

* Compute the angular momentum H, and thereby the inclination INC.

      hx(1) = x(2)*vx(3) - x(3)*vx(2)
      hx(2) = x(3)*vx(1) - x(1)*vx(3)
      hx(3) = x(1)*vx(2) - x(2)*vx(1)
      h2 = hx(1)**2 + hx(2)**2 + hx(3)**2
      h = sqrt(h2)
      if(hx(3).gt.h) then                 ! Hal's fix
        hx(3) = h
        hx(1) = 0.0d0
        hx(2) = 0.0d0
      endif
      inc = acos(hx(3)/h)

* Compute longitude of ascending node CAPOM and the argument of
* latitude u.
      fac = sqrt(hx(1)**2 + hx(2)**2)/h

      if( (fac.lt. TINY ) .or. (inc.eq.0.0d0) ) then ! Hal's fix
         capom = 0.d0
         u = atan2(x(2),x(1))
         if(abs(inc - PI).lt. 10.d0*TINY) u = -u
      else
         capom = atan2(hx(1),-hx(2))        
         u = atan2(x(3)/sin(inc),x(1)*cos(capom)+x(2)*sin(capom))
      endif

      if(capom .lt. 0.d0) capom = capom + TWOPI
      if(u .lt. 0.d0) u = u + TWOPI

*  Compute the radius R and velocity squared V2, and the dot
*  product RDOTV, the energy per unit mass ENERGY .

      r = sqrt(x(1)**2+x(2)**2+x(3)**2)
      v2 = vx(1)**2 + vx(2)**2 + vx(3)**2
      v = sqrt(v2)
      vdotr = x(1)*vx(1)+x(2)*vx(2)+x(3)*vx(3)
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
             face =(a-r)/(a*e)

c... Apr. 16/93 : watch for case where face is slightly outside unity
             if ( face .gt. 1.d0) then
                cape = 0.d0
             else
                if ( face .gt. -1.d0) then
                   cape = acos( face )
                else
                   cape = PI
                endif
             endif

            if ( vdotr .lt. 0.d0 ) cape = TWOPI - cape
          cw = (cos( cape) -e)/(1.d0 - e*cos(cape))
          sw = sqrt(1.d0 - e*e)*sin(cape)/(1.d0 - e*cos(cape))
          w = atan2(sw,cw)
          if(w .lt. 0.d0) w = w + TWOPI
        else
          e = 0.d0
          w = u
          cape = u
        endif

        capm = cape - e*sin (cape)
        omega = u - w
        if(omega .lt. 0.d0) omega = omega + TWOPI
        omega = omega - int(omega/(TWOPI))*TWOPI        

      endif
***
***
c HYPERBOLA
      if(ialpha .eq. +1) then

        a = +0.5d0*gmsum/energy  
        fac = h2/(gmsum*a)

          if (fac .gt. TINY) then
           e = sqrt ( 1.d0 + fac )
          tmpf = (a+r)/(a*e)
            if(tmpf.lt.1.0d0) then
               tmpf = 1.0d0
            endif
          capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
          if ( vdotr .lt. 0.d0 ) capf = - capf
          cw = (e - cosh(capf))/(e*cosh(capf) - 1.d0 )
          sw = sqrt(e*e - 1.d0)*sinh(capf)/(e*cosh(capf) - 1.d0 )
          w = atan2(sw,cw)
          if(w .lt. 0.d0) w = w + TWOPI
        else
c we only get here if a hyperbola is essentially a parabola
c so we calculate e and w accordingly to avoid singularities
          e = 1.d0
          tmpf = 0.5d0*h2/gmsum
          w = acos(2.d0*tmpf/r -1.d0)
          if ( vdotr .lt. 0.d0) w = TWOPI - w
          tmpf = (a+r)/(a*e)
          capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
        endif

        capm = e * sinh(capf) - capf
        omega = u - w
        if(omega .lt. 0.d0) omega = omega + TWOPI
        omega = omega - int(omega/(TWOPI))*TWOPI        
      endif
***
***
c PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)
      if(ialpha .eq. 0) then
        a =  0.5d0*h2/gmsum  
        e = 1.d0
        w = acos(2.d0*a/r -1.d0)
        if ( vdotr .lt. 0.d0) w = TWOPI - w
        tmpf = tan(0.5d0 * w)
        capm = tmpf* (1.d0 + tmpf*tmpf/3.d0)
        omega = u - w
        if(omega .lt. 0.d0) omega = omega + TWOPI
        omega = omega - int(omega/(TWOPI))*TWOPI        
      endif
***
***
      return
      end    ! orbel_xv2el
c------------------------------------------------------------------

