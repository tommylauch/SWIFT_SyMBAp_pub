c*************************************************************************
c                            HELIO_LINDRIFT_SYMBAP.F
c*************************************************************************
c This subroutine takes a linear drift due to mometum of Sun
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 vxb           ==>  velocity in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 xh            ==>  initial position in helio coord 
c                                       (real arrays)
c             Output:
c                 xh            ==>  final position in helio coord 
c                                       (real arrays)
c                 ptx           ==> momentum of sun: tp's need this   
c                                       (real scalars)
c
c Remarks: Bases on Martin's code h2.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 1/8/97

      subroutine helio_lindrift_symbap(nbod,mass,vxb,dt,xh,ptx)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt,vxb(3,nbod)

c...  Inputs and Outputs:
      real*8 xh(3,nbod)

c...  Outputs Only: 
      real*8 ptx(3)

c...  Internals:
      integer n

c----
c...  Executable code 

      ptx(:) = mass(2)*vxb(:,2)

      do n=3,nbod
         ptx(:) = ptx(:) + mass(n)*vxb(:,n)
      enddo

      ptx(:) = ptx(:)/mass(1)

      do n=2,nbod
         if(mass(n).ne.0.0d0) then
            xh(:,n) = xh(:,n) + ptx(:)*dt
         endif
      enddo

      return
      end       ! helio_lindrift_symbap
c---------------------------------------------------
