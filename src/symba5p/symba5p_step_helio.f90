!*************************************************************************
!                            SYMBA5P_STEP_HELIO.F
!*************************************************************************
! This subroutine takes a step in helio coord.  
! Does a KICK than a DRIFT than a KICK.
! ONLY DOES MASSIVE PARTICLES
!             Input:
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxh           ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 xh            ==>  final position in helio coord 
!                                       (real arrays)
!                 vxh           ==>  final velocity in helio coord 
!                                       (real arrays)
! Remarks: Based on helio_step_pl.f but does not pass the intermediate
!          positions and velocities back for the TP to use.
! Authors:  Hal Levison 
! Date:    3/20/97
! Last revision: 12/13/00

subroutine symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,        &
                              xh,vxh,dt)
use swift_mod
use coord_interface
use helio_interface
use mvs_interface
use symba5p_interface, except_this_one => symba5p_step_helio
implicit none

integer(ik), intent(in)    :: nbod,nbodm
real(rk), intent(in)       :: mass(:),dt,j2rp2,j4rp4

integer(ik), intent(inout) :: i1st

real(rk), intent(inout)    :: xh(:,:),vxh(:,:)

integer(ik)                :: i1stloc
real(rk)                   :: dth,axh(3,NTPMAX)
real(rk), save             :: vxb(3,NTPMAX)
real(rk)                   :: ptxb(3),ptxe(3)                           ! Not used here

!...  Executable code

   dth = 0.5_rk*dt
   i1stloc = i1st
   if (i1st.eq.0) then
!...     Convert vel to bery to jacobi coords
       call coord_vh2b(nbod,mass,vxh,vxb)
       i1st = 1_ik                                                      ! turn this off
   endif

!...  Do the linear drift due to momentum of the Sun
   call helio_lindrift(nbod,mass,vxb,dth,xh,ptxb)

!...  Get the accelerations in helio frame. if frist time step
   call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,xh,axh)
      i1stloc = 0_ik

!...  Apply a heliocentric kick for a half dt 
   call kickvh(nbod,vxb,axh,dth)

!...  Drift in helio coords for the full step 
   call helio_drift(nbod,mass,xh,vxb,dt)

!...  Get the accelerations in helio frame. if frist time step
   call symba5p_helio_getacch(i1stloc,nbod,nbodm,mass,j2rp2,j4rp4,xh,axh)

!...  Apply a heliocentric kick for a half dt 
   call kickvh(nbod,vxb,axh,dth)

!...  Do the linear drift due to momentum of the Sun
   call helio_lindrift(nbod,mass,vxb,dth,xh,ptxe)

!...  convert back to helio velocities
   call coord_vb2h(nbod,mass,vxb,vxh)

return
end subroutine symba5p_step_helio
