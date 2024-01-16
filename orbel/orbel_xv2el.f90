!**********************************************************************
!                      ORBEL_XV2EL.F
!**********************************************************************
!     PURPOSE:  Given the cartesian position and velocity of an orbit,
!       compute the osculating orbital elements.
!       input:
!            x,y,z    ==>  position of object (real scalars)
!            vx,vy,vz ==>  velocity of object (real scalars)
!            gmsum       ==> G*(M1+M2) (real scalar)
!       Output:
!           ialpha   ==> conic section type ( see PURPOSE, integer scalar)
!           a        ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            inc      ==> inclination  (real scalar)
!            capom    ==> longitude of ascending node (real scalar)
!           omega    ==> argument of perihelion (real scalar)
!           capm     ==> mean anomoly(real scalar)
!     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
!     REMARKS:  If the inclination INC is less than TINY, we
!       arbitrarily choose the longitude of the ascending node LGNODE
!       to be 0.0 (so the ascending node is then along the X axis).  If 
!       the  eccentricity E is less than SQRT(TINY), we arbitrarily
!       choose the argument of perihelion to be 0.
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  May 8,1992.
!     REVISIONS: 12/8/2011
!**********************************************************************

subroutine orbel_xv2el(x,vx,gmsum,ialpha,a,e,inc,capom,omega,capm)
implicit none
use swift_mod

real(rk), intent(in)     :: x(:),vx(:),gmsum

integer(ik), intent(out) :: ialpha
real(rk), intent(out)    :: a,e,inc,capom,omega,capm

real(rk)                 :: hx(3),h2,h,r,v2,v,vdotr,energy
real(rk)                 :: fac,face,cape,capf,tmpf,cw,sw,w,u

!...  Executable code 

! Compute the angular momentum H, and thereby the inclination INC.
   hx(1) = x(2)*vx(3) - x(3)*vx(2)
   hx(2) = x(3)*vx(1) - x(1)*vx(3)
   hx(3) = x(1)*vx(2) - x(2)*vx(1)
   h2 = dot_product(hx,hx)
   h = sqrt(h2)
   if (hx(3).gt.h) then                                                     ! Hal's fix
      hx(3) = h
      hx(1) = 0.0_rk
      hx(2) = 0.0_rk
   endif
   inc = acos(hx(3)/h)

! Compute longitude of ascending node CAPOM and the argument of
! latitude u.
   fac = sqrt(hx(1)**2 + hx(2)**2)/h
   if ( (fac.lt. TINY ) .or. (inc.eq.0.0_rk) ) then ! Hal's fix
      capom = 0.0_rk
      u = atan2(x(2),x(1))
      if (abs(inc - PI).lt. 10.0_rk*TINY) u = -u
   else
      capom = atan2(hx(1),-hx(2))
      u = atan2(x(3)/sin(inc),x(1)*cos(capom)+x(2)*sin(capom))
   endif
   if (capom .lt. 0.0_rk) capom = capom + TWOPI
   if (u .lt. 0.0_rk) u = u + TWOPI

!  Compute the radius R and velocity squared V2, and the dot
!  product RDOTV, the energy per unit mass ENERGY .
   r = sqrt(dot_product(x,x))
   v2 = dot_product(vx,vx)
   v = sqrt(v2)
   vdotr = dot_product(x,vx)
   energy = 0.5_rk*v2 - gmsum/r

!  Determine type of conic section and label it via IALPHA
   if (abs(energy*r/gmsum) .lt. sqrt(TINY)) then
      ialpha = 0_ik
   else
      if (energy .lt. 0.0_rk) ialpha = -1 
      if (energy .gt. 0.0_rk) ialpha = +1
   endif

! Depending on the conic type, determine the remaining elements

! ELLIPSE :
   if (ialpha.eq.-1) then
      a = -0.5_rk*gmsum/energy  
      fac = 1.0_rk-h2/(gmsum*a)
      if (fac.gt.TINY) then
         e = sqrt(fac)
         face = (a-r)/(a*e)
!... Apr. 16/93 : watch for case where face is slightly outside unity
         if (face.gt.1.0_rk) then
            cape = 0.0_rk
         else
            if (face.gt.-1.0_rk) then
               cape = acos(face)
            else
               cape = PI
            endif
         endif
         if (vdotr.lt.0.0_rk) cape = TWOPI-cape
         cw = (cos(cape)-e)/(1.0_rk-e*cos(cape))
         sw = sqrt(1.0_rk-e**2)*sin(cape)/(1.0_rk-e*cos(cape))
         w = atan2(sw,cw)
         if (w.lt.0.0_rk) w = w+TWOPI
      else
         e = 0.0_rk
         w = u
         cape = u
      endif
      capm = cape - e*sin(cape)
      omega = u - w
      if (omega .lt. 0.0_rk) omega = omega + TWOPI
      omega = omega - int(omega/(TWOPI))*TWOPI        
   endif

! HYPERBOLA
   if (ialpha.eq.1_ik) then
      a = 0.5_rk*gmsum/energy  
      fac = h2/(gmsum*a)
      if (fac .gt. TINY) then
         e = sqrt (1.0_rk+fac)
         tmpf = (a+r)/(a*e)
         if (tmpf.lt.1.0_rk) then
            tmpf = 1.0_rk
         endif
         capf = log(tmpf + sqrt(tmpf*tmpf -1.0_rk))
         if (vdotr .lt. 0.0_rk) capf = -capf
         cw = (e-cosh(capf))/(e*cosh(capf)-1.0_rk)
         sw = sqrt(e**2-1.0_rk)*sinh(capf)/(e*cosh(capf)-1.0_rk)
         w = atan2(sw,cw)
         if (w .lt. 0.0_rk) w = w+TWOPI
      else
! we only get here if a hyperbola is essentially a parabola
! so we calculate e and w accordingly to avoid singularities
         e = 1.0_rk
         tmpf = 0.5_rk*h2/gmsum
         w = acos(2.0_rk*tmpf/r-1.0_rk)
         if (vdotr .lt. 0.0_rk) w = TWOPI-w
         tmpf = (a+r)/(a*e)
         capf = log(tmpf+sqrt(tmpf*tmpf-1.0_rk))
      endif
      capm = e*sinh(capf)-capf
      omega = u-w
      if (omega .lt. 0.0_rk) omega = omega+TWOPI
      omega = omega-int(omega/(TWOPI))*TWOPI        
   endif

! PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)
   if (ialpha.eq.0_ik) then
      a =  0.5_rk*h2/gmsum  
      e = 1.0_rk
      w = acos(2.0_rk*a/r-1.0_rk)
      if (vdotr.lt.0.0_rk) w = TWOPI-w
      tmpf = tan(0.5_rk*w)
      capm = tmpf*(1.0_rk+ mpf**2/3.0_rk)
      omega = u-w
      if (omega.lt.0.0_rk) omega = omega+TWOPI
      omega = omega-int(omega/(TWOPI))*TWOPI        
   endif

return
end subroutine orbel_xv2el
