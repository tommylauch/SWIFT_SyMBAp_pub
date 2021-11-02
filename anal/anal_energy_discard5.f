c*************************************************************************
c                          ANAL_ENERGY_DISCARD5.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
c
c      Input:
c            iflg          ==>  use to turn on/off energy calculation
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            nbodm         ==>  Location of last massive body(int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh,yh,zh      ==>  current position in heliocentric coord 
c                               (real arrays)
c            vxh,vyh,vzh   ==>  current velocity in heliocentric coord 
c                               (real arrays)
c
c      Output:
c            ke            ==>  kinetic energy
c            pot           ==>  potential energy
c            energy        ==>  Total energy
c            eltot         ==>  components of total angular momentum
c                               (real array)
c
c Remarks: Based on anal_energy
c Authors:  Hal Levison
c Date:  12/16/06
c Last revision:  

      subroutine anal_energy_discard5(iflg,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,ke,pot,energy,eltot)

      include '../swift.inc'

c...  Inputs: 
      integer iflg,nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Output
      real*8 energy,eltot(3),ke,pot

c...  Internals
      logical leuse

      data leuse/.true./
      save leuse

c----
c...  Executable code

      if(iflg.lt.0) then
         leuse = .false.
         return               !  <==== NOTE
      else if(iflg.gt.0) then
         leuse = .true.
         return               !  <==== NOTE
      endif

c...  iflg = 0

      if(leuse) then
         call anal_energy_mtiny(nbod,nbodm,mass,j2rp2,j4rp4,xh,yh,zh,
     &        vxh,vyh,vzh,ke,pot,energy,eltot)
      else
         ke = 0.0d0
         pot = 0.0d0
         energy = 0.0d0
         eltot(1)=0.0d0
         eltot(2)=0.0d0
         eltot(3)=0.0d0
      endif


      return	
      end      ! anal_energy_discard5
c-----------------------------------------------------------------------

