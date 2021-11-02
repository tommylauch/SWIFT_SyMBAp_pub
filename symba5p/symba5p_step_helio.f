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
c Remarks: Based on helio_step_pl.f but does not pass the intermediate
c          positions and velocities back for the TP to use.
c Authors:  Hal Levison 
c Date:    3/20/97
c Last revision: 12/13/00

      subroutine symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(NTPMAX),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

c...  Internals:
      integer i1stloc
      real*8 dth 
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      real*8 ptxe,ptye,ptze

      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      i1stloc = i1st
      if(i1st.eq.0) then
c...      Convert vel to bery to jacobi coords
          call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)
          i1st = 1              ! turn this off
      endif

c...  Do the linear drift due to momentum of the Sun
      call helio_lindriftp(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxb,ptyb,ptzb)

c...  Get the accelerations in helio frame. if frist time step
      call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)
      i1stloc = 0

c...  Apply a heliocentric kick for a half dt 
      call kickvhp(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Drift in helio coords for the full step 
      call helio_driftp(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)

c...  Get the accelerations in helio frame. if frist time step
      call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvhp(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindriftp(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxe,ptye,ptze)

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      return
      end   ! symba5p_step_helio
c---------------------------------------------------------------------
