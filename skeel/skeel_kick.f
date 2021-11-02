c*************************************************************************
c                             SKEEL_KICK.F
c*************************************************************************
c Do a skeel kick
c
c             Input:
c                 mpl           ==>  mass of planet (real Scalar)
c                 xpl,ypl,zpl   ==>  Heliocentric position of planet 
c                                     (real Scalars)
c                 xtp,ytp,ztp   ==>  Heliocentric position of TP
c                                     (real Scalars)
c              vxtp,vytp,vztp   ==>  Heliocentric velocity of TP
c                                     (real Scalars)
c                        ri     ==>  Square of Radius of shell  (real scalar)
c                        dt     ==>  timestep  (real scalar)
c
c            Output:
c                 xtp,ytp,ztp   ==>  Heliocentric position of TP
c                                     (real Scalars)
c
c Remarks: Uses Man Hoi Lee's force
c Authors:  Hal Levison 
c Date:   9/23/96
c Last revision: 1/23/97

      subroutine skeel_kick(mpl,xpl,ypl,zpl,xtp,ytp,ztp,vxtp,
     &     vytp,vztp,ri,dt,sgn)

      include '../swift.inc'
      include 'skeel.inc'

c...  Inputs Only: 
      real*8 mpl,dt,ri,sgn
      real*8 xpl,ypl,zpl
      real*8 xtp,ytp,ztp

c...  Inputs & Outputs Only: 
      real*8 vxtp,vytp,vztp

c...  Internals: 
      real*8 ax,ay,az,r2,fac,rr,rim1,ris,r

c----
c...  Executable code 


c...  calculate the accelerations

      r2 = (xtp-xpl)**2 +  (ytp-ypl)**2 +  (ztp-zpl)**2 

      rim1 = ri*RSHELL*RSHELL

      if (r2.lt.rim1) then
         fac = 0.0d0
      else if (r2.lt.ri) then
         ris = sqrt(ri)
         r = sqrt(r2)
         rr = (ris-r)/(ris*(1.0-RSHELL))
         fac = mpl * (r2**(-1.5d0)) * 
     &        ( 1.0d0 - 3.0d0*rr*rr + 2.0d0*(rr**3))
      else
         fac = mpl * (r2**(-1.5d0))
      endif

      ax = -fac*(xtp-xpl)
      ay = -fac*(ytp-ypl)
      az = -fac*(ztp-zpl)

c...  apply the kick

      vxtp = vxtp + ax*dt*sgn
      vytp = vytp + ay*dt*sgn
      vztp = vztp + az*dt*sgn

      return
      end      ! skeel_kick.f
c--------------------------------------------------------------
