c*************************************************************************
c                        SYMBA5P_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 nbodm       ==>  Location of last massive body(int scalar)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c                 mtiny       ==>  Small mass  (real array)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks: Based on helio_getacch.f, but does not include the forces of
c          an body B on body A, if body B and A are having an encounter.
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 17/8/20

      subroutine symba5p_getacch(nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4,mtiny
      real*8 xh(nbod),yh(nbod),zh(nbod)
      integer ielst(2,NENMAX),ielc

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX) 
      integer i,j,ie,ij_tm,ij
      real*8 dx,dy,dz,rji2,faci,facj,irij3
c---
c...  Executable code 

c...  Zero things
      do i=1,nbod
         axh(i) = 0.0
         ayh(i) = 0.0
         azh(i) = 0.0
      enddo
c...  now the third terms

      ij_tm=((nbodm-1)*(nbodm-2))/2

!$OMP PARALLEL
!$OMP& REDUCTION(+:axh,ayh,azh)
!$OMP& PRIVATE(i,j,ij,ie,dx,dy,dz,rji2,irij3,faci,facj)
!$OMP& SHARED(ij_tm,nbod,nbodm,mass,xh,yh,zh,
!$OMP& ielc,ielst)
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

            axh(j) = axh(j) + (-faci*dx)
            ayh(j) = ayh(j) + (-faci*dy)
            azh(j) = azh(j) + (-faci*dz)

            axh(i) = axh(i) + facj*dx
            ayh(i) = ayh(i) + facj*dy
            azh(i) = azh(i) + facj*dz
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

         axh(j) = axh(j) + (-faci*dx)
         ayh(j) = ayh(j) + (-faci*dy)
         azh(j) = azh(j) + (-faci*dz)

         axh(i) = axh(i) + facj*dx
         ayh(i) = ayh(i) + facj*dy
         azh(i) = azh(i) + facj*dz
      enddo
!$OMP END DO NOWAIT
c...  Now subtract off anyone in an encounter
!$OMP DO
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         dx = xh(j) - xh(i)
         dy = yh(j) - yh(i)
         dz = zh(j) - zh(i)
         rji2 = dx*dx + dy*dy + dz*dz

         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3

         axh(j) = axh(j) + faci*dx
         ayh(j) = ayh(j) + faci*dy
         azh(j) = azh(j) + faci*dz

         axh(i) = axh(i) + (-facj*dx)
         ayh(i) = ayh(i) + (-facj*dy)
         azh(i) = azh(i) + (-facj*dz)
      enddo
!$OMP END DO
!$OMP END PARALLEL

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
!$OMP PARALLEL
!$OMP& PRIVATE(i)
!$OMP& SHARED(nbod,axh,ayh,azh,aoblx,aobly,aoblz)
!$OMP DO
         do i = 2,nbod
            if(mass(i).ne.0.0d0) then
               axh(i) = axh(i) + aoblx(i)
               ayh(i) = ayh(i) + aobly(i)
               azh(i) = azh(i) + aoblz(i)
            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      endif

      return
      end      ! symba5p_getacch

c---------------------------------------------------------------------
