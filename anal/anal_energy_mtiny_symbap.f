c*************************************************************************
c                          ANAL_ENERGY_MTINY_SYMBAP.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            nbodm         ==>  Location of last massive body(int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh            ==>  current position in heliocentric coord 
c                               (real array)
c            vxh           ==>  current velocity in heliocentric coord 
c                               (real array)
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

      subroutine anal_energy_mtiny_symbap(nbod,nbodm,mass,j2rp2,j4rp4,
     &           xh,vxh,ke,pot,energy,eltot)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(:),j2rp2,j4rp4
      real*8 xh(:,:),vxh(:,:)

c...  Output
      real*8 energy,eltot(:),ke,pot

c...  Internals
      real*8 elx(3)
      real*8 xx(3),rr2,oblpot,msys,irh(NTPMAX),ir3h(NTPMAX)
      real*8 xb(3,NTPMAX)       ! Used NTPMAX for symba
      real*8 vxb(3,NTPMAX)
      integer i,j

c----
c...  Executable code 

      call coord_h2b_symbap(nbod,mass,xh,vxh,xb,vxb,msys)   

      eltot(1) = (xb(2,nbod)*vxb(3,nbod)-xb(3,nbod)*vxb(2,nbod))
      eltot(2) = (xb(3,nbod)*vxb(1,nbod)-xb(1,nbod)*vxb(3,nbod))
      eltot(3) = (xb(1,nbod)*vxb(2,nbod)-xb(2,nbod)*vxb(1,nbod))
      eltot = eltot*mass(nbod)

      ke = 0.5*mass(nbod)*(vxb(1,nbod)**2+vxb(2,nbod)**2+vxb(3,nbod)**2)
      pot= 0.d0

      do i=1,nbodm
         elx(1) = xb(2,i)*vxb(3,i)-xb(3,i)*vxb(2,i)
         elx(2) = xb(3,i)*vxb(1,i)-xb(1,i)*vxb(3,i)
         elx(3) = xb(1,i)*vxb(2,i)-xb(2,i)*vxb(1,i)
         elx = elx*mass(i)
         eltot(:) = eltot(:) + elx(:)
         
         ke = ke + 0.5*mass(i)*(vxb(1,i)**2 + vxb(2,i)**2 + vxb(3,i)**2)
         do j = i+1,nbod
            xx(:) = xb(:,i) - xb(:,j)
            rr2 = xx(1)**2 + xx(2)**2 + xx(3)**2
            if((mass(i).ne.0.0d0).and.(mass(j).ne.0.0d0)) then
               pot = pot - mass(i)*mass(j)/(sqrt(rr2))
            endif
         enddo
      enddo

      do i=nbodm+1,nbod-1
         ke = ke + 0.5*mass(i)*(vxb(1,i)**2 + vxb(2,i)**2 + vxb(3,i)**2)
         elx(1) = xb(2,i)*vxb(3,i)-xb(3,i)*vxb(2,i)
         elx(2) = xb(3,i)*vxb(1,i)-xb(1,i)*vxb(3,i)
         elx(3) = xb(1,i)*vxb(2,i)-xb(2,i)*vxb(1,i)
         elx = elx*mass(i)
         eltot(:) = eltot(:) + elx(:)
      enddo

      if(j2rp2.ne.0.0d0) then
         call getacch_ir3_symbap(nbod,2,xh,ir3h,irh)
         call obl_pot_symbap(nbod,mass,j2rp2,j4rp4,xh,irh,oblpot)
         pot = pot + oblpot
      endif

      energy = ke + pot

      return	
      end      ! anal_energy_symbap
c-----------------------------------------------------------------------

