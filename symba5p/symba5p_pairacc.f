c*************************************************************************
c                        SYMBA5P_PAIRACC.F
c*************************************************************************
c Mutual acceleration of all pairs (i<=nbodm, j>i), heliocentric positions
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  location of the last massive body (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  position in helio coord (real array)
c             Output:
c                 acc           ==>  acceleration, acc(:,1) = 0 (real array)
c
c Remarks: SoA copies so the inner j loop vectorizes; reduction arrays are
c          sized nbod (not NTPMAX).
c Authors:  TCH Lau
c Date:    11/10/26

      subroutine symba5p_pairacc(nbod,nbodm,mass,xh,acc)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod,nbodm
      real*8 mass(*),xh(3,*)

c...  Outputs:
      real*8 acc(3,*)

c...  Internals:
      integer i,j
      real*8 x(nbod),y(nbod),z(nbod),ax(nbod),ay(nbod),az(nbod)
      real*8 xi,yi,zi,mi,axi,ayi,azi,dx,dy,dz,r2,irij3,faci,facj

c----
c...  Executable code

      do j=1,nbod
         x(j) = xh(1,j)
         y(j) = xh(2,j)
         z(j) = xh(3,j)
         ax(j) = 0.0d0
         ay(j) = 0.0d0
         az(j) = 0.0d0
      enddo

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1)
!$OMP& PRIVATE(i,j,xi,yi,zi,mi,axi,ayi,azi)
!$OMP& PRIVATE(dx,dy,dz,r2,irij3,faci,facj)
!$OMP& SHARED(nbod,nbodm,mass,x,y,z)
!$OMP& REDUCTION(+:ax,ay,az)
      do i=2,nbodm
         xi = x(i)
         yi = y(i)
         zi = z(i)
         mi = mass(i)
         axi = 0.0d0
         ayi = 0.0d0
         azi = 0.0d0
!$OMP SIMD REDUCTION(+:axi,ayi,azi)
!$OMP& PRIVATE(dx,dy,dz,r2,irij3,faci,facj)
         do j=i+1,nbod
            dx = x(j) - xi
            dy = y(j) - yi
            dz = z(j) - zi
            r2 = dx*dx + dy*dy + dz*dz
            irij3 = 1.0d0/(r2*sqrt(r2))
            faci = mi*irij3
            facj = mass(j)*irij3
            ax(j) = ax(j) - faci*dx
            ay(j) = ay(j) - faci*dy
            az(j) = az(j) - faci*dz
            axi = axi + facj*dx
            ayi = ayi + facj*dy
            azi = azi + facj*dz
         enddo
         ax(i) = ax(i) + axi
         ay(i) = ay(i) + ayi
         az(i) = az(i) + azi
      enddo
!$OMP END PARALLEL DO

      do j=1,nbod
         acc(1,j) = ax(j)
         acc(2,j) = ay(j)
         acc(3,j) = az(j)
      enddo
      acc(:,1) = 0.0d0

      return
      end      ! symba5p_pairacc
c--------------------------------------------------------------------------
