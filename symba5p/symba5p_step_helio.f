c*************************************************************************
c                            SYMBA5P_STEP_HELIO.F
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
c                 xh            ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh           ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh            ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh           ==>  final velocity in helio coord 
c                                       (real arrays)
c Remarks: Based on helio_step_pl.f but does not pass the intermediate
c          positions and velocities back for the TP to use.
c Authors:  Hal Levison 
c Date:    3/20/97
c Last revision: 12/13/00

      subroutine symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,vxh,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(3,nbod)
      real*8 vxh(3,nbod)

c...  Internals:
      integer i1stloc
      real*8 dth,axh(3,NTPMAX),vxb(3,NTPMAX),msys
      real*8 ptxb(3)            ! Not used here
      real*8 ptxe(3)
      save vxb     ! Note this !!
c----
c...  Executable code 

      dth = 0.5d0*dt

      i1stloc = i1st
      if(i1st.eq.0) then
c...      Convert vel to bery to jacobi coords
          call coord_vh2b_symbap(nbod,mass,vxh,vxb,msys)
          i1st = 1              ! turn this off
      endif

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift_symbap(nbod,mass,vxb,dth,xh,ptxb)

c...  Get the accelerations in helio frame. if frist time step
      call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,axh)
      i1stloc = 0

c...  Apply a heliocentric kick for a half dt 
      call kickvh_symbap(nbod,vxb,axh,dth)

c..   Drift in helio coords for the full step 
      call helio_drift_symbap(nbod,mass,xh,vxb,dt)

c...  Get the accelerations in helio frame. if frist time step
      call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,axh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_symbap(nbod,vxb,axh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift_symbap(nbod,mass,vxb,dth,xh,ptxe)

c...  convert back to helio velocities
      call coord_vb2h_symbap(nbod,mass,vxb,vxh)

      return
      end   ! symba5p_step_helio
c---------------------------------------------------------------------
