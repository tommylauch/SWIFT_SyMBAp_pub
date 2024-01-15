c*************************************************************************
c                            DISCARD_MASS_REORDER5_SYMBAP.F
c*************************************************************************
c Remove a massive body
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ip            ==>  planets to remove (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>   position in helio coord 
c                                    (real arrays)
c                 vxh           ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 rhill         ==>  size of a planet's hill's sphere.
c                                    (real array)
c                 isperih       ==> heliocentric peri flags. (real array)
c             Output:
c                 ip            ==>  planets to remove (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>   position in helio coord 
c                                    (real arrays)
c                 vxh           ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 rhill         ==>  size of a planet's hill's sphere.
c                                    (real array)
c                 isperih       ==> heliocentric peri flags. (real array)
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    1/2/97
c Last revision: 5/13/99

      subroutine discard_mass_reorder5_symbap(ip,nbod,mass,xh,vxh,
     &                                        rpl,rhill,isperih)

      include '../swift.inc'

c...  Inputs: 
      integer ip

c...  Input and Output
      integer nbod
      real*8 mass(nbod),xh(3,nbod)
      real*8 vxh(3,nbod),rpl(nbod)
      real*8 rhill(nbod)
      integer isperih(nbod)

c...  internal
      integer i,j

c-----
c...  Executable code 

      do i=ip,nbod-1
         xh(:,i) = xh(:,i+1)
         vxh(:,i) = vxh(:,i+1)
         mass(i) = mass(i+1)
         rpl(i) = rpl(i+1)
         rhill(i) = rhill(i+1)
         isperih(i) = isperih(i+1)
      enddo
      nbod = nbod - 1

      return
      end
