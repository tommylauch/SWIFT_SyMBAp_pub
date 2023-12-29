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
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(3,nbod)

c...  Outputs:
      real*8 axh(3,nbod)
                
c...  Internals:
      integer i,j
      real*8 aoblx(3,NTPMAX),axhl(3,NTPMAX)
      real*8 ir3h(NTPMAX),irh(NTPMAX)
      real*8 dx(3),rji2,faci,facj,irij3

      save axhl     ! Note this !!

c----
c...  Executable code 

      if(iflg.eq.0) then
         axhl = 0.0
c...     now the third terms
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& REDUCTION(+:axhl)
!$OMP& PRIVATE(i,j,dx,rji2,irij3,faci,facj)
!$OMP& SHARED(nbod,nbodm,mass,xh)
!$OMP DO COLLAPSE(2)
      do i=2,nbodm
         do j=i+1,nbod
            dx(:) = xh(:,j) - xh(:,i)
            rji2 = dx(1)**2 + dx(2)**2 + dx(3)**2

            irij3 = 1.0d0/(rji2*sqrt(rji2))
            faci = mass(i)*irij3
            facj = mass(j)*irij3

            axhl(:,j) = axhl(:,j) - faci*dx(:)
            axhl(:,i) = axhl(:,i) + facj*dx(:)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      endif
c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3_symbap(nbod,2,xh,ir3h,irh)
         call obl_acc_symbap(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
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
