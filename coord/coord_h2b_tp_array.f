c***********************************************************************
c	                    COORD_H2B_TP_ARRAY.F
c***********************************************************************
*     PURPOSE: Converts test part from Heliocentric to Barycentric coords.
*     ARGUMENTS:  Input is 
*                              ntp ==> number of test part (<= NTPMAX)
*                                              (integer)
*		        xht,yht,zht ==> heliocentric particle coords
*                                          (real array)
*		     vxht,vyht,vzht ==> heliocentric particle velocities
*                                             (real array)
*		        xhs,yhs,zhs ==> bary coords of the Sun
*                                          (real scalar)
*		     vxhs,vyhs,vzhs ==> bary vel of the Sun
*                                          (real scalar)
*                 Returned are
*                       xbt,ybt,zbt ==> bary. particle positions
*                                            (real array)
*                    vxbt,vybt,vzbt ==> bary. particle velocities
*                                            (real array)
*       
*     Authors:  Hal Levison
*     ALGORITHM: Obvious 
*     WRITTEN:  2/18/92
*     REVISIONS:

      subroutine coord_h2b_tp_array(ntp,xht,vxht,xs,vxs,xbt,vxbt)

      include '../swift.inc'

c...  Inputs: 
      integer ntp
      real*8 xht(3,NTPMAX),vxht(3,NTPMAX)
      real*8 xs(3),vxs(3)

c...  Outputs:
      real*8 xbt(3,NTPMAX),vxbt(3,NTPMAX)

c...  Internals:
      integer i

c----
c...  Executable code 
      do i=1,ntp
         xbt(:,i) = xht(:,i) + xs
         vxbt(:,i) = vxht(:,i) + vxs
      enddo

      return
      end     ! coord_h2b_tp_array
c--------------------------------------------------------------------------

