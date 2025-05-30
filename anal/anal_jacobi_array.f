c*************************************************************************
c                          ANAL_JACOBI_array.F
c*************************************************************************
c Calculates the Jacobi integral of a particle given the position and velocity
c and ang. velocity vector of the Sun and Planet
c
c      Input:
c            msun          ==>   Mass of the Sun (real scalar)
c            mpl           ==>   Mass of the Planet (real scalar)
c            omegax,..y,z  ==>   vector components of planet-sun orbit's
c                                 angular velocity vector(real scalars)
c            xb,yb,zb      ==>   Barycentric position of particle 
c                                   (real scalars)
c            vxb,vyb,vzb   ==>   Barycentric vel of particle (real scalars)
c            xsb,ysb,zsb   ==>   Barycentric position of Sun (real scalars)
c            xplb,yplb,zplb ==>  Barycentric position of planet (real scalars)
c
c       Output:
c            jacobi         ==>  Value of the jacobi constant
c
c Remarks: 
c Authors:  Martin Duncan
c Date:   April 13/93 
c Last revision:  3/4/93 HFL 

      subroutine anal_jacobi_array(msun,mpl,omegax,xb,
     &      vxb,xsb,xplb,jacobi)

      include '../swift.inc'

c...  Inputs: 
      real*8 msun,mpl,omegax(3),xb(3),vxb(3)
      real*8 xsb(3),xplb(3)

c...  Outputs:
      real*8 jacobi

c...  Internals
      real*8 rr,jx(3)

c----
c...  Executable code 

      jacobi = 0.5d0*sum(vxb**2)

      rr = sqrt(sum((xb-xsb)**2))
      jacobi = jacobi - msun/rr

      rr = sqrt(sum((xb-xplb)**2))
      jacobi = jacobi - mpl/rr

      jx(1) = xb(2)*vxb(3) - xb(3)*vxb(2)
      jx(2) = xb(3)*vxb(1) - xb(1)*vxb(3)
      jx(3) = xb(1)*vxb(2) - xb(2)*vxb(1)

      jacobi = jacobi - sum(omegax*jx)

      return
      end       ! anal_jacobi_array
c-------------------------------------------------------------------------





