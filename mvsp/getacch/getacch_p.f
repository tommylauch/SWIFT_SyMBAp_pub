c*************************************************************************
c                        GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/18/93

      subroutine getacch_p(nbod,mass,j2rp2,j4rp4,xj,xh,axh)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(NPLMAX),xj(3,NPLMAX),j2rp2,j4rp4
      real*8 xh(3,NPLMAX)

c...  Outputs:
      real*8 axh(3,NPLMAX)
                
c...  Internals:
      integer i
      real*8 ir3h(NPLMAX),ir3j(NPLMAX)
      real*8 irh(NPLMAX),irj(NPLMAX)
      real*8 axh1(3,NPLMAX),axh2(3,NPLMAX),axh3(3,NPLMAX),axh0(3)
      real*8 aoblx(3,NPLMAX)

c----
c...  Executable code 

c...  get thr r^-3's
      call getacch_ir3_p(nbod,2,xj,ir3j,irj)
      call getacch_ir3_p(nbod,2,xh,ir3h,irh)

c...  calc the ah0's:  recall that they are the same for all particles
      call getacch_ah0_p(3,nbod,mass,xh,ir3h,axh0) 

c...  now the first terms
      call getacch_ah1_p(nbod,mass,xh,xj,ir3h,ir3j,axh1)

c...  now the second terms
      call getacch_ah2_array(nbod,mass,xj,ir3j,axh2)

c...  now the third terms
      call getacch_ah3_p(nbod,mass,xh,axh3)

c...  add them all together
      axh(:,1) = 0.d0
      do i=2,nbod
        axh(:,i) = axh0 + axh1(:,i) + axh2(:,i) + axh3(:,i)
      enddo

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call obl_acc_array(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
         do i = 2,nbod
            axh(:,i) = axh(:,i) + aoblx(:,i) - aoblx(:,1)
         enddo
      endif

      return
      end      ! getacch_p

c---------------------------------------------------------------------




