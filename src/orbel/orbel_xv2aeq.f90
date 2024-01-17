!***********************************************************************
!                      ORBEL_XV2AEQ.F
!***********************************************************************
!       PURPOSE:  Given the cartesian position and velocity of an orbit,
!         compute the osculating orbital elements a, e, and q only.
!       input:
!            x        ==>  position of object (real scalars)
!            vx       ==>  velocity of object (real scalars)
!            gmsum    ==> G*(M1+M2) (real scalar)
!       Output:
!           ialpha    ==> conic section type ( see PURPOSE, integer scalar)
!           a         ==> semi-major axis or pericentric distance if a parabola
!                          (real scalar)
!            e        ==> eccentricity (real scalar)
!            q        ==> perihelion distance (real scalar); q = a(1 - e)
!       ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech." 
!       REMARKS: Based on M. Duncan's orbel_xv2el.f
!        This routine is generally applied to study (hyperbolic) close 
!       encounters of test particles with planets.
!       AUTHOR:  L. Dones
!       DATE WRITTEN:  February 24, 1994
!       REVISIONS: 

subroutine orbel_xv2aeq(x,vx,gmsum,ialpha,a,e,q)
use swift_mod
implicit none

real(rk), intent(in)     :: x(:),vx(:),gmsum

integer(ik), intent(out) :: ialpha
real(rk), intent(out)    :: a,e,q

real(rk)                 :: hx(3),h2,r,v2,energy,fac

!...  Executable code

! Compute the angular momentum H, and thereby the inclination INC.
   hx(1) = x(2)*vx(3) - x(3)*vx(2)
   hx(2) = x(3)*vx(1) - x(1)*vx(3)
   hx(3) = x(1)*vx(2) - x(2)*vx(1)

   h2 = dot_product(hx,hx)

!    Compute the radius R and velocity squared V2, and the dot
!    product RDOTV, the energy per unit mass ENERGY .
   r = sqrt(dot_product(x(1:3),x(1:3)))
   v2 = dot_product(vx(1:3),vx(1:3))
   energy = 0.5_rk*v2 - gmsum/r

!    Determine type of conic section and label it via IALPHA
   if (abs(energy*r/gmsum) .lt. sqrt(TINY)) then
      ialpha = 0_ik
   else
      if (energy .lt. 0.0_rk) ialpha = -1_ik
      if (energy .gt. 0.0_rk) ialpha = +1_ik
   endif

! Depending on the conic type, determine the remaining elements

! ELLIPSE :
   if (ialpha.eq.-1) then
      a = -0.5_rk*gmsum/energy  
      fac = 1.0_rk - h2/(gmsum*a)
      if (fac.gt.TINY) then
         e = sqrt (fac)
      else
         e = 0.0_rk
      endif
      q = a*(1.0_rk-e)
   endif

! HYPERBOLA
   if (ialpha.eq.1) then
      a = 0.5_rk*gmsum/energy  
      fac = h2/(gmsum*a)
      if (fac.gt.TINY) then
         e = sqrt(1.0_rk+fac)
         q = -a*(1.0_rk-e)
!     have to insert minus sign in expression for q because this code
!      takes a > 0, even for a hyperbola
      else
! we only get here if a hyperbola is essentially a parabola
! so we calculate e accordingly to avoid singularities
         e = 1.0_rk
         q = 0.5_rk*h2/gmsum
      endif
   endif

! PARABOLA : ( NOTE - in this case "a", which is formally infinite,
!         is arbitrarily set equal to the pericentric distance q).
   if (ialpha.eq.0) then
     a =  0.5_rk*h2/gmsum  
     e = 1.0_rk
     q = a
   endif

return
end subroutine orbel_xv2aeq
