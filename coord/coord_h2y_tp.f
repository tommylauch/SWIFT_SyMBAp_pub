c***********************************************************************
c	                    COORD_H2Y_TP.F
c***********************************************************************
*     PURPOSE: Converts test part from Heliocentric to Yosemite coords.
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*                       xht,yht,zht ==> helio. particle positions
*                                            (real array)
*                    vxht,vyht,vzht ==> helio. particle velocities
*                                            (real array)
*		        xyo,yyo,zyo ==> offset yose coords 
*                                          (real scalar)
*		     vxyo,vyyo,vzyo ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*		        xyt,yyt,zyt ==> yose particle coords
*                                          (real array)
*		     vxyt,vyyt,vzyt ==> yose particle velocities
*                                             (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  9/13/02
*     REVISIONS:

	subroutine coord_h2y_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &      xyo,yyo,zyo,vxyo,vyyo,vzyo,
     &      xyt,yyt,zyt,vxyt,vyyt,vzyt)


      include '../swift.inc'

c...  Inputs: 
	integer ntp
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 xyo,yyo,zyo,vxyo,vyyo,vzyo

c...  Outputs:
	real*8 xyt(NTPMAX),yyt(NTPMAX),zyt(NTPMAX)
	real*8 vxyt(NTPMAX),vyyt(NTPMAX),vzyt(NTPMAX)

c...  Internals:
	integer i

c----
c...  Executable code 
	do i=1,ntp
	  xyt(i) = xht(i) - xyo
	  yyt(i) = yht(i) - yyo
	  zyt(i) = zht(i) - zyo
	  vxyt(i) = vxht(i) - vxyo
	  vyyt(i) = vyht(i) - vyyo
	  vzyt(i) = vzht(i) - vzyo
	enddo

	return
	end     ! coord_h2y_tp
c--------------------------------------------------------------------------

