c*************************************************************************
c                        SYMBA5P_HELIO_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 iflg        ==>  =0 calculate forces (int scalor)
c                                  =1 don't
c                 nbod        ==>  number of massive bodies (int scalor)
c                 nbodm       ==>  The last massive particle
c                                  (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh          ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh         ==>  acceleration in helio coord (real arrays)
c
c Remarks Based on helio_getacch.f
c Author:  Hal Levison  
c Date:    9/12/99
c Last revision: 11/08/13 

      subroutine symba5p_helio_getacch(iflg,nbod,nbodm,mass,
     &     j2rp2,j4rp4,xh,axh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,nbodm,iflg
      real*8 mass(*),j2rp2,j4rp4
      real*8 xh(3,*)

c...  Outputs:
      real*8 axh(3,*)
                
c...  Internals:
      integer i
      real*8 aoblx(3,NTPMAX),axhl(3,NTPMAX)
      real*8 ir3h(NTPMAX),irh(NTPMAX)

      save axhl     ! Note this !!

c----
c...  Executable code 

      if(iflg.eq.0) then
         call symba5p_pairacc(nbod,nbodm,mass,xh,axhl)
      endif
c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3_p(nbod,2,xh,ir3h,irh)
         call obl_acc_array(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
         do i=2,nbod
            axh(:,i) = axhl(:,i) + aoblx(:,i)
         enddo
      else
         do i=2,nbod
            axh(:,i) = axhl(:,i)
         enddo
      endif

      return
      end      ! symba5p_helio_getacch

c---------------------------------------------------------------------
