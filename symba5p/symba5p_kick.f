c*************************************************************************
c                             SYMBA5P_KICK.F
c*************************************************************************
c Do a symba5 kick
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 irec          ==>  recursion level  (integer scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 rhill         ==>  Hill sphere of planet (real Scalar)
c                 xh            ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb           ==>  initial velocity in bari coord 
c                                    (real arrays)
c                dt             ==>  timestep  (real scalar)
c                sgn            ==>  add or subtract the force (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c            Output:
c                 vxb           ==>  final velocity in bari coord 
c                                    (real arrays)
c
c Remarks: Uses Man Hoi's force
c Authors:  Hal Levison 
c Date:   3/20/97
c Last revision: 3/3/10

      subroutine symba5p_kick(nbod,mass,irec,iecnt,ielev,
     &                        rhill,xh,vxb,dt,sgn,ielc,ielst)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer nbod,irec
      real*8 mass(*),dt,rhill(*),sgn
      integer iecnt(*),ielev(*)
      real*8 xh(3,*)
      integer ielst(2,*),ielc

c...  Inputs and Outputs:
      real*8 vxb(3,*)

c...  Internals: 
      real*8 dx(3),fac,ris,r
      real*8 ri,rr,r2,faci,facj,ir3,rim1
      integer i,j,irm1,irecl,ie

c----
c...  Executable code 

      irm1 = irec - 1
      if(sgn.lt.0.0d0) then
         irecl = irec - 1
      else
         irecl = irec
      endif

c...  calculate the accelerations
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)

         if((ielev(i).ge.irm1) .and. (ielev(j).ge.irm1) ) then

            ri = (rhill(i)+rhill(j))**2 * RHSCALE*RHSCALE * 
     &           (RSHELL**(2*irecl))
            rim1 = ri*RSHELL*RSHELL
            
            dx(:) = xh(:,j) - xh(:,i)
            r2 = dx(1)**2 + dx(2)**2 + dx(3)**2
            ir3 = 1.0d0/(r2*sqrt(r2))

            if (r2.lt.rim1) then
               fac = 0.0d0
            else if (r2.lt.ri) then
               ris = sqrt(ri)
               r = sqrt(r2)
               rr = (ris-r)/(ris*(1.0-RSHELL))
               fac = ir3 * (1.d0 - (rr**4)*(35.d0 -
     &                      rr*(84.d0 - rr*(70.d0 - 20.d0*rr))))
            else
               fac = ir3
            endif

            if( (iecnt(i).ne.0) .and. (ielev(i).ge.irm1) ) then
               facj = mass(j)*fac
               vxb(:,i) = vxb(:,i) + facj*dx(:)*dt*sgn
            endif

            if( (iecnt(j).ne.0) .and. (ielev(j).ge.irm1) ) then
               faci = mass(i)*fac
               vxb(:,j) = vxb(:,j) - faci*dx(:)*dt*sgn
            endif

         endif

      enddo
      return
      end      ! symba5p_kick.f
c--------------------------------------------------------------