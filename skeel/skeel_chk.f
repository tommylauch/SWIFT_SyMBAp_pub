c*************************************************************************
c                            SKEEL_CHK.F
c*************************************************************************
c This subroutine checks to see if there are encounters
c
c             Input:
c                 mpl           ==>  mass of the planet (real scalar)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real scalar)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real scalar)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real scalar)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real scalar)
c                 dt            ==>  time step  (real scalor)
c                 r2crit        ==>  critical distence (real scalor)
c             Output:
c                 icflg         ==> ecounters? = 1 Yes
c                                              =  0 No (integer scalar)  
c
c Remarks: Based on RMVS3_CHK.F
c Authors:  Hal Levison 
c Date:    9/24/96
c Last revision: 

      subroutine skeel_chk(mpl,xh,yh,zh,vxh,vyh,vzh,xht,yht,
     &       zht,vxht,vyht,vzht,dt,r2crit,icflg)

      include '../swift.inc'
      include 'skeel.inc'

c...  Inputs: 
      real*8 mpl,xh,yh,zh,dt
      real*8 xht,yht,zht
      real*8 vxh,vyh,vzh
      real*8 vxht,vyht,vzht
      real*8 r2crit

c...  Outputs
       integer icflg

c...  Internals
       integer iflag
       real*8 xr,yr,zr,vxr,vyr,vzr
       real*8 r2critp

c-----
c...  Executable code 

       xr = xht - xh
       yr = yht - yh
       zr = zht - zh
       vxr = vxht - vxh
       vyr = vyht - vyh
       vzr = vzht - vzh

       r2critp = 0.0            ! dummy so we can use rmvs routine

       call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &      r2crit,r2critp,iflag)

       if(iflag.gt.0) then
          icflg = 1
       else
          icflg = 0
       endif
		 
       return
       end                      ! skeel_chk
c------------------------------------------------------

