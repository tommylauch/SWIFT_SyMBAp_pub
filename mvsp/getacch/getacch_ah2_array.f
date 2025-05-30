c*************************************************************************
c                        GETACCH_AH2.F
c*************************************************************************
c This subroutine calculates the 2nd term of acceleration 
c on the massive particles in the HELIOCENTRIC frame. Only in one direction. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 xj,yj,zj    ==>  position in jacobi coord (real array)
c                 ir3j        ==> inv radii in jacobi coord (real array)
c             Output:
c                 axh2,ayh2,azh2 ==>  2nd term acceleration in helio coord 
c                                    (real array)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/3/93

      subroutine getacch_ah2_array(nbod,mass,xj,ir3j,axh2)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),ir3j(nbod)
      real*8 xj(3,nbod)

c...  Outputs:
      real*8 axh2(3,nbod)
                
c...  Internals:
      integer i
      real*8 etaj, fac

c----
c...  Executable code 

      axh2(:,1) = 0.d0
      axh2(:,2) = 0.d0

      etaj = mass(1)
      do i=3,nbod
         etaj = etaj + mass(i-1)
         fac = mass(i)*mass(1)*ir3j(i)/etaj
         axh2(:,i) = axh2(:,i-1) + fac*xj(:,i)
      enddo

      return
      end     ! getacch_ah2
c---------------------------------------------------------------------




