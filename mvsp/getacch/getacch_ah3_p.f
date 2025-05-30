c*************************************************************************
c                        GETACCH_AH3_P.F
c*************************************************************************
c This subroutine calculates the 3rd term acceleration on the massive particles
c in the HELIOCENTRIC frame. This term is the direct cross terms
c             Input:
c                 nbod          ==>  number of massive bodies (int scalor)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  position in heliocentric coord 
c                                   (real arrays)
c             Output:
c                 axh3,ayh3,azh3 ==>  3rd term of acceleration in helio coord 
c                                     (real arrays)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 11/21/96

      subroutine getacch_ah3_p(nbod,mass,xh,axh3)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),xh(3,nbod)

c...  Outputs:
      real*8 axh3(3,nbod)

c...  Internals:
      integer i,j
      real*8 dx(3),rji2,faci,facj,irij3

c------
c...  Executable code

      axh3 = 0.d0

      do i=2,nbod-1
         do j=i+1,nbod
             dx = xh(:,j)-xh(:,i)
             rji2 = sum(dx**2)

             irij3 = 1.0d0/(rji2*sqrt(rji2))
             faci = mass(i)*irij3
             facj = mass(j)*irij3

             axh3(:,j) = axh3(:,j) - faci*dx
             axh3(:,i) = axh3(:,i) + facj*dx
         enddo
      enddo

      return
      end     ! getacch_ah3_p

