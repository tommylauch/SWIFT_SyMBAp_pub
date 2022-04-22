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
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks Based on helio_getacch.f
c Author:  Hal Levison  
c Date:    9/12/99
c Last revision: 11/08/13 

      subroutine symba5p_helio_getacch(iflg,nbod,nbodm,mass,
     &     j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,nbodm,iflg
      real*8 mass(NTPMAX),j2rp2,j4rp4
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)

c...  Outputs:
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
                
c...  Internals:
      integer i,j,ij_tm,ij
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 axhl(NTPMAX),ayhl(NTPMAX),azhl(NTPMAX)
      real*8 ir3h(NTPMAX),irh(NTPMAX)
      real*8 dx,dy,dz,rji2,faci,facj,irij3

      save axhl,ayhl,azhl     ! Note this !!

c----
c...  Executable code 

      if(iflg.eq.0) then
         do i=1,nbod
            axhl(i) = 0.0
            ayhl(i) = 0.0
            azhl(i) = 0.0
         enddo

c...     now the third terms
      ij_tm=((nbodm-1)*(nbodm-2))/2

!$OMP PARALLEL
!$OMP& REDUCTION(+:axhl,ayhl,azhl)
!$OMP& PRIVATE(i,j,ij,dx,dy,dz,rji2,irij3,faci,facj)
!$OMP& SHARED(ij_tm,nbod,nbodm,mass,xh,yh,zh)
!$OMP DO COLLAPSE(2)
      do i=2,nbodm
         do j=nbodm+1,nbod
            dx = xh(j) - xh(i)
            dy = yh(j) - yh(i)
            dz = zh(j) - zh(i)
            rji2 = dx*dx + dy*dy + dz*dz

            irij3 = 1.0d0/(rji2*sqrt(rji2))
            faci = mass(i)*irij3
            facj = mass(j)*irij3

            axhl(j) = axhl(j) + (-faci*dx)
            ayhl(j) = ayhl(j) + (-faci*dy)
            azhl(j) = azhl(j) + (-faci*dz)

            axhl(i) = axhl(i) + facj*dx
            ayhl(i) = ayhl(i) + facj*dy
            azhl(i) = azhl(i) + facj*dz
         enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO
      do ij=0,(ij_tm-1)
         i=ij/(nbodm-2)+2 !i goes from 2
         j=mod(ij,(nbodm-2))+3  !j goes from 3
         if(j.le.i) then
            i = nbodm-i+2 !i goes to nbodm-1
            j = nbodm-j+3 !j goes to nbodm
         endif

         dx = xh(j) - xh(i)
         dy = yh(j) - yh(i)
         dz = zh(j) - zh(i)
         rji2 = dx*dx + dy*dy + dz*dz

         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3

         axhl(j) = axhl(j) + (-faci*dx)
         ayhl(j) = ayhl(j) + (-faci*dy)
         azhl(j) = azhl(j) + (-faci*dz)

         axhl(i) = axhl(i) + facj*dx
         ayhl(i) = ayhl(i) + facj*dy
         azhl(i) = azhl(i) + facj*dz
      enddo
!$OMP END DO
!$OMP END PARALLEL

      endif
c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
!$OMP PARALLEL
!$OMP& PRIVATE(i)
!$OMP& SHARED(nbod,axh,ayh,azh,axhl,ayhl,azhl,
!$OMP& aoblx,aobly,aoblz)
!$OMP DO
         do i = 2,nbod
            axh(i) = axhl(i) + aoblx(i)
            ayh(i) = ayhl(i) + aobly(i)
            azh(i) = azhl(i) + aoblz(i)
         enddo
!$OMP END DO
!$OMP END PARALLEL
      else
!$OMP PARALLEL
!$OMP& PRIVATE(i)
!$OMP& SHARED(nbod,axh,ayh,azh,axhl,ayhl,azhl)
!$OMP DO
         do i = 2,nbod
            axh(i) = axhl(i)
            ayh(i) = ayhl(i)
            azh(i) = azhl(i)
         enddo
!$OMP END DO
!$OMP END PARALLEL
      endif

      return
      end      ! symba5p_helio_getacch

c---------------------------------------------------------------------
