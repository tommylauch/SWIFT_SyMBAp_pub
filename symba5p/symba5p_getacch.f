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
c                 xh          ==>  position in heliocentric coord (real arrays)
c                 mtiny       ==>  Small mass  (real array)
c                ielc         ==>  number of encounters (integer scalar)
c                ielst        ==>  list of ecnounters (2D integer array)
c             Output
c                 axh         ==>  acceleration in helio coord (real arrays)
c
c Remarks: Based on helio_getacch.f, but does not include the forces of
c          an bodxy B on body A, if body B and A are having an encounter.
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 17/8/20

      subroutine symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,axh,
     &                           mtiny,ielc,ielst)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4,mtiny
      real*8 xh(3,nbod)
      integer ielst(2,NENMAX),ielc

c...  Outputs:
      real*8 axh(3,nbod)

c...  Internals:
      real*8 aoblx(3,NTPMAX)
      real*8 ir3h(NTPMAX),irh(NTPMAX) 
      integer i,j,ie
      real*8 dx(3),rji2,faci,facj,irij3
c---
c...  Executable code 

c...  Zero things
      axh = 0.d0
c...  now the third terms
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& REDUCTION(+:axh)
!$OMP& PRIVATE(i,j,ie,dx,rji2,irij3,faci,facj)
!$OMP& SHARED(nbod,nbodm,mass,xh,ielc,ielst)
!$OMP DO COLLAPSE(2)
      do i=2,nbodm
         do j=i+1,nbod
            dx(:) = xh(:,j) - xh(:,i)
            rji2 = dx(1)**2 + dx(2)**2 + dx(3)**2

            irij3 = 1.0d0/(rji2*sqrt(rji2))
            faci = mass(i)*irij3
            facj = mass(j)*irij3

            axh(:,j) = axh(:,j) - faci*dx(:)
            axh(:,i) = axh(:,i) + facj*dx(:)
         enddo
      enddo
!$OMP END DO NOWAIT
c...  Now subtract off anyone in an encounter
!$OMP DO
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         dx(:) = xh(:,j) - xh(:,i)
         rji2 = dx(1)**2 + dx(2)**2 + dx(3)**2

         irij3 = 1.0d0/(rji2*sqrt(rji2))
         faci = mass(i)*irij3
         facj = mass(j)*irij3

         axh(:,j) = axh(:,j) + faci*dx(:)
         axh(:,i) = axh(:,i) - facj*dx(:)
      enddo
!$OMP END DO
!$OMP END PARALLEL

c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3_symbap(nbod,2,xh,ir3h,irh)
         call obl_acc_symbap(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
         do i=2,nbod
            if(mass(i).ne.0.0d0) then
               axh(:,i) = axh(:,i) + aoblx(:,i)
            endif
         enddo
      endif

      return
      end      ! symba5p_getacch

c---------------------------------------------------------------------
