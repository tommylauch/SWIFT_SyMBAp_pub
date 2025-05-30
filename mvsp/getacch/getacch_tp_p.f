c*************************************************************************
c                        GETACCH_TP_P.F
c*************************************************************************
c This subroutine calculates the acceleration on the test particles
c in the HELIOCENTRIC frame. 
c             Input:
c                  nbod        ==>  number of massive bodies (int scalor)
c                  ntp         ==>  number of tp bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                  j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                  xh,yh,zh    ==>  massive part position in helio coord 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position in heliocentric coord 
c                                     (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c               axht,ayht,azht ==>  tp acceleration in helio coord 
c                                   (real arrays)
c
c Author:  Hal Levison  
c Date:    2/18/93
c Last revision: 2/24/94

      subroutine getacch_tp_p(nbod,ntp,mass,j2rp2,j4rp4,xh,
     &     xht,istat,axht)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX)
      real*8 mass(NPLMAX),xh(3,NPLMAX)
      real*8 xht(3,NTPMAX),j2rp2,j4rp4

c...  Outputs:
      real*8 axht(3,NTPMAX)
                
c...  Internals:
      integer i
      real*8 ir3h(NPLMAX),irh(NPLMAX)
      real*8 ir3ht(NTPMAX),irht(NTPMAX)
      real*8 axh3(3,NTPMAX)
      real*8 axh0(3)
      real*8 aoblx(3,NPLMAX),aoblxt(3,NTPMAX)

c----
c...  Executable code 


c...  get thr r^-3's  for the planets and test paricles
      call getacch_ir3_p(nbod,2,xh,ir3h,irh)
      call getacch_ir3_p(ntp,1,xht,ir3ht,irht)

c...  calc the ah0's:  recall that they are the same for all particles
      call getacch_ah0_p(2,nbod,mass,xh,ir3h,axh0) 

c...  the first terms are 0

c...  the second terms are 0

c...  now the third terms
      call getacch_ah3_tp_p(nbod,ntp,mass,xh,xht,istat,axh3)

c...  add them all together
      do i=1,ntp
         if(istat(i).eq.0) then
            axht(:,i) = axh0 + axh3(:,i)
         endif
      enddo

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call obl_acc_array(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
         call obl_acc_tp_array(ntp,istat,mass(1),j2rp2,j4rp4,xht,
     &        irht,aoblxt)
         do i = 1,ntp
            axht(:,i) = axht(:,i) + aoblxt(:,i) - aoblx(:,1)
         enddo
      endif

      return
      end      ! getacch_tp_p

c---------------------------------------------------------------------




