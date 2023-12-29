c*************************************************************************
c                            UTIL_HILLS_SYMBAP.F
c*************************************************************************
c This subroutine calculates the hill's sphere for the planets
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  initial position in helio coord 
c                                    (real array)
c                 vxh           ==>  initial velocity in helio coord 
c                                    (real array)
c             Output:
c                  r2hill       ==>  the SQUARE of the planet's hill's sphere 
c                                    (real array)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 1/6/97

      subroutine util_hills_symbap(nbod,mass,xh,vxh,r2hill) 

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      real*8 mass(nbod),xh(3,nbod)
      real*8 vxh(3,nbod)

c...  Outputs
      real*8 r2hill(nbod)

c...  Internals
      integer i
      real*8 mu,energy,ap,rhil,r,v2

c-----
c...  Executable code 

      do i=2,nbod
         if(mass(i).ne.0.0d0) then
            mu = mass(1)*mass(i)/(mass(1)+mass(i))
            r = sqrt(xh(1,i)**2+xh(2,i)**2+xh(3,i)**2)
            v2 = vxh(1,i)**2+vxh(2,i)**2+vxh(3,i)**2
            energy =-1.0d0*mass(1)*mass(i)/r + 0.5*mu*v2
            ap = -1.0d0*mass(1)*mass(i)/(2.0d0*energy)
            rhil = ap * (((mu/mass(1))/3.0)**(0.3333333333))
            r2hill(i) = rhil*rhil
         else
            r2hill(i) = 0.0d0
         endif
      enddo
      
      r2hill(1) = 0.0
      
      return
      end                       ! util_hills_symbap

c---------------------------------------------------
