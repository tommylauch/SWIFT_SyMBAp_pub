c*************************************************************************
c                            STEP_KDK_PL_P.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/24/94

      subroutine step_kdk_pl_p(i1st,nbod,mass,j2rp2,j4rp4,
     &     xh,vxh,dt)

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(3,nbod),vxh(3,nbod)

c...  Internals:
      real*8 dth 
      real*8 axh(3,NPLMAX),xj(3,NPLMAX),vxj(3,NPLMAX)

      save axh,xj     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then
c...      Convert to jacobi coords
          call coord_h2j_array(nbod,mass,xh,vxh,xj,vxj)
c...     Get the accelerations in helio frame. if frist time step
         call getacch_p(nbod,mass,j2rp2,j4rp4,xj,xh,axh)
         i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_p(nbod,vxh,axh,dth)

c...  Convert the helio. vels. to Jac. vels. (the Jac. positions are unchanged)
      call coord_vh2vj_array(nbod,mass,vxh,vxj)

c..   Drift in Jacobi coords for the full step 
      call drift_p(nbod,mass,xj,vxj,dt)

c...  After drift, compute helio. xh and vh for acceleration calculations
      call coord_j2h_array(nbod,mass,xj,vxj,xh,vxh)

c...  Get the accelerations in helio frame.
      call getacch_p(nbod,mass,j2rp2,j4rp4,xj,xh,axh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_p(nbod,vxh,axh,dth)

      return
      end   ! step_kdk_pl_p
c---------------------------------------------------------------------

