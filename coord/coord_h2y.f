c***********************************************************************
c	                    COORD_H2Y.F
c***********************************************************************
*     PURPOSE: Converts from Helio to Yosemite coords.
*     ARGUMENTS:  Input is 
*                    nbod ==> number of bodies (must be less than NBMAX)
*                             (integer)
*	             mass(*) ==>  masses (real array)
*		     xh(*),yh(*),zh(*) ==> Helio particle coords
*                                          (real array)
*		     vxh(*),vyh(*),vzh(*) ==> Helio particle velocities
*                                             (real array)
*                 Returned are
*                    xy(*),yy(*),zy(*) ==>  Yose particle positions
*                                          (real array)
*                    vxy(*),vyy(*),vzy(*) ==> Yose particle velocities
*                                            (real array)
*                    xyo,yyo,zyo,vxyo,vyyo,vzyo ==>  Yose Offset vectors 
*                                                    (real scalors)
*                    mu        ==>  `reduced' Yose mass (real scalor)
*       
*     ALGORITHM: Obvious 
*     REMARKS:  
*
*     Authors:  Hal Levison
*     WRITTEN:  9/13/02

      subroutine coord_h2y(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &     xy,yy,zy,vxy,vyy,vzy,xyo,yyo,zyo,vxyo,vyyo,vzyo,mu)
      
      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(NPLMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

c...  Outputs:
      real*8 xy(NPLMAX),yy(NPLMAX),zy(NPLMAX)
      real*8 vxy(NPLMAX),vyy(NPLMAX),vzy(NPLMAX)
      real*8 xyo,yyo,zyo,vxyo,vyyo,vzyo,mu

c...  Internals:
      integer i
      real*8 const,p1

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
         xyo = xyo + mass(i)*xh(i)
         yyo = yyo + mass(i)*yh(i)
         zyo = zyo + mass(i)*zh(i)
         vxyo = vxyo + mass(i)*vxh(i)
         vyyo = vyyo + mass(i)*vyh(i)
         vzyo = vzyo + mass(i)*vzh(i)
      enddo

      mu = mu/mass(1)
      p1 = sqrt(1.0d0 + mu)
      const = (p1 - 1.0d0)/(mu*p1)

      xyo = xyo * const / mass(1)
      yyo = yyo * const / mass(1)
      zyo = zyo * const / mass(1)
      vxyo = vxyo * const / mass(1)
      vyyo = vyyo * const / mass(1)
      vzyo = vzyo * const / mass(1)

      xy(1) = - xyo
      yy(1) = - yyo
      zy(1) = - zyo
      vxy(1) = - vxyo
      vyy(1) = - vyyo
      vzy(1) = - vzyo
      do i=2,nbod
         xy(i) = xh(i) - xyo
         yy(i) = yh(i) - yyo
         zy(i) = zh(i) - zyo
         vxy(i) = vxh(i) - vxyo
         vyy(i) = vyh(i) - vyyo
         vzy(i) = vzh(i) - vzyo
      enddo

      return
      end                       ! coord_h2y

c--------------------------------------------------------------------------

