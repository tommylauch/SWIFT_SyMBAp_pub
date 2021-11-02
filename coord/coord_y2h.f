c***********************************************************************
c	                    COORD_Y2H.F
c***********************************************************************
*     PURPOSE: Converts from Yosemite to Helio coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*                    xy(*),yy(*),zy(*) ==>  Yose particle positions
*                                          (real array)
*                    vxy(*),vyy(*),vzy(*) ==> Yose particle velocities
*                                            (real array)
*                 Returned are
*		     xh(*),yh(*),zh(*) ==> Helio particle coords
*                                          (real array)
*		     vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                             (real array)
*                    xyo,yyo,zyo,vxyo,vyyo,vzyo ==>  Yose Offset vectors (real scalors)
*       
*     ALGORITHM: Obvious 
*     REMARKS:  
*
*     Authors:  Hal Levison
*     WRITTEN:  9/13/02

      subroutine coord_y2h(nbod,mass,xy,yy,zy,vxy,vyy,vzy,
     &     xh,yh,zh,vxh,vyh,vzh,xyo,yyo,zyo,vxyo,vyyo,vzyo)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(NPLMAX)
      real*8 xy(NPLMAX),yy(NPLMAX),zy(NPLMAX)
      real*8 vxy(NPLMAX),vyy(NPLMAX),vzy(NPLMAX)

c...  Outputs:
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
      real*8 xyo,yyo,zyo,vxyo,vyyo,vzyo

c...  Internals:
      integer i
      real*8 const,mu,p1

c---- 
c...  Executable code 

      mu = 0.0d0
      xyo = 0.0d0
      yyo = 0.0d0
      zyo = 0.0d0
      vxyo = 0.0d0
      vyyo = 0.0d0
      vzyo = 0.0d0
      do i=2,nbod
         mu = mu + mass(i)
         xyo = xyo + mass(i)*xy(i)
         yyo = yyo + mass(i)*yy(i)
         zyo = zyo + mass(i)*zy(i)
         vxyo = vxyo + mass(i)*vxy(i)
         vyyo = vyyo + mass(i)*vyy(i)
         vzyo = vzyo + mass(i)*vzy(i)
      enddo

      mu = mu/mass(1)
      p1 = sqrt(1.0d0 + mu)
      const = (p1 - 1.0d0) / mu

      xyo = xyo * const / mass(1)
      yyo = yyo * const / mass(1)
      zyo = zyo * const / mass(1)
      vxyo = vxyo * const / mass(1)
      vyyo = vyyo * const / mass(1)
      vzyo = vzyo * const / mass(1)

      xh(1) = 0.0d0
      yh(1) = 0.0d0
      zh(1) = 0.0d0
      vxh(1) = 0.0d0
      vyh(1) = 0.0d0
      vzh(1) = 0.0d0
      do i=2,nbod
         xh(i) = xy(i) + xyo
         yh(i) = yy(i) + yyo
         zh(i) = zy(i) + zyo
         vxh(i) = vxy(i) + vxyo
         vyh(i) = vyy(i) + vyyo
         vzh(i) = vzy(i) + vzyo
      enddo

      return
      end                       ! coord_y2h

c--------------------------------------------------------------------------

