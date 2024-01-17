c*************************************************************************
c                        DRIFT_DAN_SYMBAP.F
c*************************************************************************
c This subroutine does the Danby and decides which vbles to use
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 x0            ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx0           ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 dt0           ==>  time step
c             Output:
c                 x0            ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx0           ==>  final position in jacobi coord 
c                                       (real scalars)
c                 iflg          ==>  integer flag (zero if satisfactory)
c                                    (non-zero if nonconvergence)
c
c Authors:  Hal Levison & Martin Duncan  
c Date:    2/10/93
c Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged

      subroutine drift_dan_symbap(mu,x0,vx0,dt0,iflg)

      include '../../swift.inc'

c...  Inputs Only: 
      real*8 mu,dt0

c...  Inputs and Outputs:
      real*8 x0(*),vx0(*)

c...  Output
      integer iflg

c...  Internals:
      real*8 x(3),vx(3),dt
      real*8 f,g,fdot,c1,c2
      real*8 c3,gdot
      real*8 u,alpha,fp,r0,v0s
      real*8 a,asq,en
      real*8 dm,ec,es,esq,xkep
      real*8 fchk,s,c

c----
c...  Executable code 

c...  Set dt = dt0 to be sure timestep is not altered while solving
c...  for new coords.
      dt = dt0
      iflg = 0
      r0 = sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)
      v0s = vx0(1)**2 + vx0(2)**2 + vx0(3)**2
      u = x0(1)*vx0(1) + x0(2)*vx0(2) + x0(3)*vx0(3)
      alpha = 2.0*mu/r0 - v0s

      if (alpha.gt.0.d0) then
         a = mu/alpha
         asq = a*a
         en = sqrt(mu/(a*asq))
         ec = 1.0d0 - r0/a
         es = u/(en*asq)
         esq = ec*ec + es*es
         dm = dt*en - int(dt*en/TWOPI)*TWOPI
         dt = dm/en
         if((dm*dm .gt. 0.16d0) .or. (esq.gt.0.36d0)) goto 100

         if(esq*dm*dm .lt. 0.0016) then
            call drift_kepmd(dm,es,ec,xkep,s,c)
            fchk = (xkep - ec*s +es*(1.-c) - dm)

            if(fchk*fchk .gt. DANBYB) then
               iflg = 1
               return
            endif

            fp = 1. - ec*c + es*s
            f = (a/r0) * (c-1.) + 1.
            g = dt + (s-xkep)/en
            fdot = - (a/(r0*fp))*en*s
            gdot = (c-1.)/fp + 1.

            x(1:3) = x0(1:3)*f + vx0(1:3)*g
            vx(1:3) = x0(1:3)*fdot + vx0(1:3)*gdot

            x0(1:3) = x(1:3)
            vx0(1:3) = vx(1:3)
            iflg = 0
            return
         endif
      endif
             
100   call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

      if(iflg .eq.0) then
           f = 1.0 - (mu/r0)*c2
           g = dt - mu*c3
           fdot = -(mu/(fp*r0))*c1
           gdot = 1. - (mu/fp)*c2

           x(1:3) = x0(1:3)*f + vx0(1:3)*g
           vx(1:3) = x0(1:3)*fdot + vx0(1:3)*gdot

           x0(1:3) = x(1:3)
           vx0(1:3) = vx(1:3)
      endif

      return
      end   ! drift_dan_symbap
