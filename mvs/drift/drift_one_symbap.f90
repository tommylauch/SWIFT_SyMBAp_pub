c*************************************************************************
c                        DRIFT_ONE_SYMBAP.F
c*************************************************************************
c This subroutine does the danby-type drift for one particle, using 
c appropriate vbles and redoing a drift if the accuracy is too poor 
c (as flagged by the integer iflg).
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mu            ==>  mass of central body (real scalar) 
c                 x,y,z         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx,vy,vz      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dt            ==>  time step
c             Output:
c                 x,y,z         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx,vy,vz      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 iflg          ==>  integer (zero for successful step)
c
c Authors:  Hal Levison & Martin Duncan 
c Date:    2/10/93
c Last revision: 2/10/93
c

      subroutine drift_one_symbap(mu,x,vx,dt,iflg)

      include '../../swift.inc'

c...  Inputs Only: 
      real*8 mu,dt

c...  Inputs and Outputs:
      real*8 x(3),vx(3)

c...  Output
      integer iflg
      
c...  Internals:
      integer i
      real*8 dttmp

c----
c...  Executable code 
      call drift_dan_symbap(mu,x,vx,dt,iflg)
      if(iflg .ne. 0) then
         do i = 1,10
            dttmp = dt/10.d0
            call drift_dan_symbap(mu,x,vx,dttmp,iflg)
            if(iflg .ne. 0) return
         enddo
      endif
      return
      end    ! drift_one_symbap
c-------------------------------------------------------------------
