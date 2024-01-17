c*************************************************************************
c                            SYMBA5P_CHK.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 rhill         ==>  Radius of hill sphere (real array)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ip1,ip2       ==>  The two bodies to check (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step  (real scalor)
c                 irec          ==>  current recursion level (int scalar)
c             Output:
c                 icflg         ==> ecounter?  = 1 Yes
c                                              =  0 No (integer scalar)  
c                 svdotr        ==> = .true. if i,j are receding
c                                   = .false is approaching
c                                     (logical*1 scalar)
c
c Remarks: Based on plh_chk.f.  Same as symba_chk.f
c Authors:  Hal Levison
c Date:   3/20/97
c Last revision: 
c

      subroutine symba5p_chk(rhill,nbod,ip1,ip2,mass,xh,
     &     vxh,dt,irec,icflg,svdotr)

      include '../swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs: 
      integer nbod,irec,ip1,ip2
      real*8 mass(*),xh(3,*),dt
      real*8 vxh(3,*),rhill(*)

c...  Outputs
      integer icflg
      logical*1 svdotr

c...  Internals
      real*8 r2crit,r2critp,rcrit
      real*8 xr(3),vxr(3)
      real*8 vdotr

c-----
c...  Executable code 

      rcrit = (rhill(ip1)+rhill(ip2)) * RHSCALE * (RSHELL**(irec))
      r2crit = rcrit*rcrit
      r2critp = -1.0d0          ! not used here
      
      xr(:) = xh(:,ip2) - xh(:,ip1)
      vxr(:) = vxh(:,ip2) - vxh(:,ip1)
      call rmvs_chk_ind_symbap(xr,vxr,dt,r2crit,r2critp,icflg)
      
      vdotr = xr(1)*vxr(1) + xr(2)*vxr(2) + xr(3)*vxr(3)
      svdotr = (vdotr.lt.0.0d0)

      return
      end                       ! symba5p_chk
c------------------------------------------------------


