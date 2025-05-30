c*************************************************************************
c                        GETACCH_A3_TP_P.F
c*************************************************************************
c This subroutine calculates the 3rd term acceleration on the test particles
c in the HELIOCENTRIC frame. This term is the direct cross terms
c             Input:
c                 nbod          ==>  number of massive bodies (int scalor)
c                 ntp           ==>  number of test particles (int scalor)
c                 mass          ==>  mass of massive bodies (real array)
c                 xh,yh,zh      ==>  massive position in heliocentric coord 
c                                   (real arrays)
c                 xht,yht,zht   ==>  tp position in heliocentric coord 
c                                   (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c                 axh3,ayh3,azh3 ==>  3rd term of acceleration in helio coord 
c                                     (real arrays)
c
c Author:  Hal Levison  
c Date:    2/18/93
c Last revision: 

      subroutine getacch_ah3_tp_p(nbod,ntp,mass,xh,xht,istat,axh3)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(3,NPLMAX)
      real*8 xht(3,NTPMAX)

c...  Outputs:
      real*8 axh3(3,NTPMAX)

c...  Internals:
      integer i,j
      real*8 dx(3),rji2,fac,irij3

c------
c...  Executable code

      axh3 = 0.d0

      do j=1,ntp
         if(istat(j).eq.0) then
            do i=2,nbod
               dx = xht(:,j)-xh(:,i)
               rji2 = sum(dx**2)

               irij3 = 1.0d0/(rji2*sqrt(rji2))
               fac = mass(i)*irij3

               axh3(:,j) = axh3(:,j) - fac*dx
            enddo
         endif
      enddo

      return
      end     ! getacch_ah3_tp_p
c--------------------------------------------------------------------------

