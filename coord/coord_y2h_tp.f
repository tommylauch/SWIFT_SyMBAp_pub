c***********************************************************************
c	                    COORD_Y2H_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Yosemite to Heliocentric coords.
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*		        xyt,yyt,zyt ==> yose particle coords
*                                          (real array)
*		     vxyt,vyyt,vzyt ==> yose particle velocities
*                                             (real array)
*		        xyo,yyo,zyo ==> offset yose coords 
*                                          (real scalar)
*		     vxyo,vyyo,vzyo ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*                       xht,yht,zht ==> helio. particle positions
*                                            (real array)
*                    vxht,vyht,vzht ==> helio. particle velocities
*                                            (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  9/13/02
*     REVISIONS:

      subroutine coord_y2h_tp(ntp,xyt,yyt,zyt,vxyt,vyyt,vzyt,
     &     xyo,yyo,zyo,vxyo,vyyo,vzyo,
     &     xht,yht,zht,vxht,vyht,vzht)

      include '../swift.inc'

c...  Inputs: 
      integer ntp
      real*8 xyt(NTPMAX),yyt(NTPMAX),zyt(NTPMAX)
      real*8 vxyt(NTPMAX),vyyt(NTPMAX),vzyt(NTPMAX)
      real*8 xyo,yyo,zyo,vxyo,vyyo,vzyo

c...  Outputs:
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

c...  Internals:
      integer i

c---- 
c...  Executable code 
      do i=1,ntp
         xht(i) = xyt(i) + xyo
         yht(i) = yyt(i) + yyo
         zht(i) = zyt(i) + zyo
         vxht(i) = vxyt(i) + vxyo
         vyht(i) = vyyt(i) + vyyo
         vzht(i) = vzyt(i) + vzyo
      enddo

      return
      end                       ! coord_y2h_tp
c--------------------------------------------------------------------------

