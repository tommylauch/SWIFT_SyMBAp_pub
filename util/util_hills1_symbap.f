c*************************************************************************
c                            UTIL_HILLS1_symbap.F
c*************************************************************************
c This subroutine calculates the hill's sphere for the planets
c
c             Input:
c                 msun          ==>  mass of sun (real scalar)
c                 mpl           ==>  mass of sun (real scalar)
c                 xh            ==>  position of pl in helio coord 
c                                    (real scalars)
c                 vxh           ==>  velocity of pl in helio coord 
c                                    (real scalars)
c             Output:
c                  rhill        ==>  the radius of planet's hill's sphere 
c                                    (real scalar)
c
c
c Remarks: Based on util_hill
c Authors:  Hal Levison 
c Date:    1/8/97
c Last revision: 

      subroutine util_hills1_symbap(msun,mpl,xh,vxh,rhill) 

      include '../swift.inc'

c...  Inputs: 
      real*8 msun,mpl,xh(3),vxh(3)

c...  Outputs
      real*8 rhill

c...  Internals
      real*8 mu,energy,ap,r,v2

c-----
c...  Executable code 

      mu = msun*mpl/(msun+mpl)
      r = sqrt( xh(1)**2 + xh(2)**2 + xh(3)**2 )
      v2 =  vxh(1)**2 + vxh(2)**2 + vxh(3)**2
      energy = -1.0d0*msun*mpl/r + 0.5*mu*v2
      ap = -1.0d0*msun*mpl/(2.0d0*energy)
      rhill = ap * (((mu/msun)/3.0)**(0.3333333333))
      
      return
      end                       ! util_hills1_symbap

c---------------------------------------------------
