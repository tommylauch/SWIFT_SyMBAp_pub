c*************************************************************************
c                            IO_DISCARD_MERGE_SYMBAP
c*************************************************************************
c Write out information about a merger.
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 ip1,ip2       ==>  planets to merge (real scalar)
c                 m1            ==>  Mass of pl 1 (real scalar)
c                 r1            ==>  Radius of pl 1 (real scalar)
c                 x1            ==>  current position of pl 1 in helio coord 
c                                    (real arrays)
c                 vx1           ==>  current velocity of pl 1 in helio coord 
c                                    (real arrays)
c                 m2            ==>  Mass of pl 2 (real scalar)
c                 r2            ==>  Radius of pl 2 (real scalar)
c                 x2            ==>  current position of pl 2 in helio coord 
c                                    (real arrays)
c                 vx2           ==>  current velocity of pl 2 in helio coord 
c                                    (real arrays)
c                 mn            ==>  Mass of new pl  (real scalar)
c                 rn            ==>  Radius of new pl (real scalar)
c                 xn            ==>  current position of new pl in helio coord 
c                                    (real arrays)
c                 vxn           ==>  current velocity of new pl in helio coord 
c                                    (real arrays)
c                 nleft         ==>  number of active test bodies(int scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 

      subroutine io_discard_merge_symbap(time,ip1,ip2,m,r,x,vx,
     &                            mn,rn,xn,vxn)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer ip1,ip2
      real*8 time
      real*8 m(*),r(*),x(3,*),vx(3,*)
      real*8 mn,rn,xn(*),vxn(*)

c...  Internals
      integer ierr,iu

c----
c...  Executable code 

      iu = 40

      call io_open(iu,'discard_mass.out','append','FORMATTED',ierr)

      write(iu,1000) time
 1000 format(1x,1p1e23.16,'  2')

      write(iu,2000) ip1,m(1),r(1)
 2000 format('-1',1x,i7,1x,2(1p1e23.16,1x))
      write(iu,3000) x(1,1),x(2,1),x(3,1)
 3000 format(3(1p1e23.16,1x))
      write(iu,3000) vx(1,1),vx(2,1),vx(3,1)

      write(iu,2000) ip2,m(2),r(2)
      write(iu,3000) x(1,2),x(2,2),x(3,2)
      write(iu,3000) vx(1,2),vx(2,2),vx(3,2)

      write(iu,4000) ip1,mn,rn
 4000 format('+1',1x,i7,1x,2(1p1e23.16,1x))
      write(iu,3000) xn(1),xn(2),xn(3)
      write(iu,3000) vxn(1),vxn(2),vxn(3)

      close(unit = iu)
      return
      end                       ! io_discard_merge_symbap.f
c--------------------------------------------------------------------------

