!*************************************************************************
!                            SYMBA5P_MERGE.F
!*************************************************************************
! This subroutine checks to see if there are encounters
!             Input:
!                 t             ==>  current time (real scalar)
!                 dt            ==>  time step (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 ip1,ip2       ==>  The two bodies to check (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb           ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 ireci         ==>  current recursion level (int scalar)
!                 irecl         ==>  maximum recursion level (int scalar)
!                 svdotrold     ==>  old radial velocity test
!                                   = .true. if i,j are receding
!                                   = .false is approaching
!                                     (logical*1 scalar)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!             rhill             ==>  Hill sphere of planet (real Scalar)
!             eoff              ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
!             Output:  Changed only if a Megrer happens
!                 mass          ==>  mass of bodies (real array)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb           ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!             rhill             ==>  Hill sphere of planet (real Scalar)
!             eoff              ==>  Energy offset (real scalar)
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
! Remarks: 
! Authors:  Hal Levison
! Date:   1/2/97
! Last revision: 12/16/09

subroutine symba5p_merge(t,dt,ip1,ip2,mass,xh,vxb,svdotrold,           &
                         rpl,mergelst,mergecnt,rhill,eoff,ielc,ielst)
use swift_mod
use symba5p_mod
use discard_interface
use orbel_interface
use util_interface
implicit none

integer(ik), intent(in)    :: ip1,ip2
real(rk), intent(in)       :: t,dt
logical(ik), intent(in)    :: svdotrold

real(rk), intent(inout)    :: mass(:),xh(:,:),vxb(:,:),eoff
real(rk), intent(inout)    :: rpl(:),rhill(:)
integer(ik), intent(inout) :: ielst(:,:),ielc
integer(ik), intent(inout) :: mergelst(:,:),mergecnt

integer(ik)                :: ip1l,ip2l,ialpha
real(rk)                   :: xr(3),vxr(3),vdotr,tcross2
real(rk)                   :: rlim,rlim2,rr2,massc,a,e,peri,dt2

!...  Executable code 

   xr(:) = xh(:,ip2) - xh(:,ip1)
   rr2 = dot_product(xr,xr)
   rlim = rpl(ip1)+rpl(ip2)

   if (rlim.eq.0.0_rk) return  ! <======  NOTE !!!!!

   rlim2 = rlim**2

   if (rlim2.ge.rr2) then
      ip1l = ip1
      ip2l = ip2 
!$OMP CRITICAL (MERGE)
      call discard_mass_merge5p_mtiny(t,ip1l,ip2l,mass,xh,vxb,rpl,     &
                                      eoff,ielc,ielst)
      mergecnt = mergecnt+1
      mergelst(1,mergecnt) = ip1l
      mergelst(2,mergecnt) = ip2l
      rhill(ip2l) = 0.0_rk
      call util_hills1(mass(1),mass(ip1l),xh(1:3,ip1l),                &
                       vxb(1:3,ip1l),rhill(ip1l))
!$OMP END CRITICAL (MERGE)
      return      !   <=== NOTE !!!!!!!!!
   endif

   vxr(:) = vxb(:,ip2) - vxb(:,ip1)
   vdotr = dot_product(xr,vxr)

   if ( svdotrold .and. (vdotr.gt.0.0_rk) ) then
      tcross2 = rr2/(dot_product(vxr,vxr))
      dt2 = dt**2

      if (tcross2.le.dt2) then
         massc = mass(ip1)+mass(ip2)
         call orbel_xv2aeq(xr,vxr,massc,ialpha,a,e,peri)
         if (peri.lt.rlim) then
            ip1l = ip1
            ip2l = ip2 
!$OMP CRITICAL (MERGE)
            call discard_mass_merge5p_mtiny(t,ip1l,ip2l,mass,xh,vxb,rpl,&
                                            eoff,ielc,ielst)
            mergecnt = mergecnt + 1
            mergelst(1,mergecnt) = ip1l
            mergelst(2,mergecnt) = ip2l
            rhill(ip2l) = 0.0_rk
            call util_hills1(mass(1),mass(ip1l),xh(1:3,ip1l),       &
                             vxb(1:3,ip1l),rhill(ip1l))
!$OMP END CRITICAL (MERGE)
         endif
      endif
   endif

return
end subroutine symba5p_merge
