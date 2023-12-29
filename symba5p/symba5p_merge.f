c*************************************************************************
c                            SYMBA5P_MERGE.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 t             ==>  current time (real scalar)
c                 dt            ==>  time step (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 ip1,ip2       ==>  The two bodies to check (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb           ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 ireci         ==>  current recursion level (int scalar)
c                 irecl         ==>  maximum recursion level (int scalar)
c                 svdotrold     ==>  old radial velocity test
c                                   = .true. if i,j are receding
c                                   = .false is approaching
c                                     (logical*1 scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c             mergelst          ==>  list of mergers (int array)
c             mergecnt          ==>  count of mergers (int array)
c             rhill             ==>  Hill sphere of planet (real Scalar)
c             eoff              ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c
c             Output:  Changed only if a Megrer happens
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb           ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c             mergelst          ==>  list of mergers (int array)
c             mergecnt          ==>  count of mergers (int array)
c             rhill             ==>  Hill sphere of planet (real Scalar)
c             eoff              ==>  Energy offset (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c
c Remarks: 
c Authors:  Hal Levison
c Date:   1/2/97
c Last revision: 12/16/09
c

      subroutine symba5p_merge(t,dt,nbod,nbodm,ip1,ip2,mass,xh,vxb,
     &     ireci,irecl,svdotrold,iecnt,rpl,mergelst,mergecnt,rhill,
     &     eoff,ielc,ielst)


      include '../swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs: 
      integer nbod,nbodm,ireci,irecl,ip1,ip2
      real*8 t,dt
      logical*1 svdotrold

c...  Inputs and Outputs:
      real*8 mass(NTPMAX),xh(3,NTPMAX),vxb(3,NTPMAX),eoff
      real*8 rpl(NTPMAX),rhill(NTPMAX)
      integer iecnt(NTPMAX)
      integer mergelst(2,NTPMAX),mergecnt
      integer ip1l,ip2l
      integer ielst(2,NENMAX),ielc

c...  Outputs

c...  Internals
      integer ialpha
      real*8 xr(3),vxr(3),vdotr,tcross2
      real*8 rlim,rlim2,rr2,massc,a,e,peri,dt2

c-----
c...  Executable code 

      xr(1) = xh(1,ip2) - xh(1,ip1)
      xr(2) = xh(2,ip2) - xh(2,ip1)
      xr(3) = xh(3,ip2) - xh(3,ip1)
      rr2 = xr(1)**2 + xr(2)**2 + xr(3)**2

      rlim = rpl(ip1)+rpl(ip2)

      if(rlim.eq.0.0d0) RETURN  ! <======  NOTE !!!!!

      rlim2 = rlim*rlim

      if(rlim2.ge.rr2) then
         ip1l = ip1
         ip2l = ip2 
!$OMP CRITICAL (MERGE)
         call discard_mass_merge5p_mtiny(t,nbod,nbodm,ip1l,ip2l,
     &                  mass,xh,vxb,rpl,eoff,ielc,ielst,NENMAX)
         mergecnt = mergecnt + 1
         mergelst(1,mergecnt) = ip1l
         mergelst(2,mergecnt) = ip2l
         rhill(ip2l) = 0.0d0
         call util_hills1_symbap(mass(1),mass(ip1l),xh(:,ip1l),
     &                           vxb(:,ip1l),rhill(ip1l))
!$OMP END CRITICAL (MERGE)
         return      !   <=== NOTE !!!!!!!!!
      endif

      vxr(1) = vxb(1,ip2) - vxb(1,ip1)
      vxr(2) = vxb(2,ip2) - vxb(2,ip1)
      vxr(3) = vxb(3,ip2) - vxb(3,ip1)
      vdotr = xr(1)*vxr(1) + xr(2)*vxr(2) + xr(3)*vxr(3)

      if( svdotrold .and. (vdotr.gt.0.0d0)) then

         tcross2 = rr2/(vxr(1)**2+vxr(2)**2+vxr(3)**2)
         dt2 = dt*dt

         if(tcross2.le.dt2) then
            massc = mass(ip1) + mass(ip2)
            call orbel_xv2aeq_symbap(xr,vxr,massc,ialpha,a,e,peri)
            if( peri.lt.rlim) then
               ip1l = ip1
               ip2l = ip2 
!$OMP CRITICAL (MERGE)
               call discard_mass_merge5p_mtiny(t,nbod,nbodm,ip1l,ip2l,
     &                        mass,xh,vxb,rpl,eoff,ielc,ielst,NENMAX)
               mergecnt = mergecnt + 1
               mergelst(1,mergecnt) = ip1l
               mergelst(2,mergecnt) = ip2l
               rhill(ip2l) = 0.0d0
               call util_hills1_symbap(mass(1),mass(ip1l),xh(:,ip1l),
     &                                 vxb(:,ip1l),rhill(ip1l))
!$OMP END CRITICAL (MERGE)
            endif
         endif
      endif

      return
      end                       ! symba5p_merge
c------------------------------------------------------

