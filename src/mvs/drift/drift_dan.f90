!*************************************************************************
!                         DRIFT_DAN_SYMBAP.F
!*************************************************************************
! This subroutine does the Danby and decides which vbles to use
!              Input:
!                  nbod          ==>  number of massive bodies (int scalar)
!                  mass          ==>  mass of bodies (real array)
!                  x0            ==>  initial position in jacobi coord 
!                                     (real scalar)
!                  vx0           ==>  initial position in jacobi coord 
!                                     (real scalar)
!                  dt0           ==>  time step
!              Output:
!                  x0            ==>  final position in jacobi coord 
!                                        (real scalars)
!                  vx0           ==>  final position in jacobi coord 
!                                        (real scalars)
!                  iflg          ==>  integer flag (zero if satisfactory)
!                                     (non-zero if nonconvergence)
! Authors:  Hal Levison & Martin Duncan  
! Date:    2/10/93
! Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged

subroutine drift_dan(mu,x0,vx0,dt0,iflg)
use swift_mod
use mvs_interface, except_this_one => drift_dan
implicit none

real(rk), intent(in)     :: mu,dt0

real(rk), intent(inout)  :: x0(:),vx0(:)

integer(ik), intent(out) :: iflg

real(rk)                :: x(3),vx(3),dt
real(rk)                :: f,g,fdot,c1,c2
real(rk)                :: c3,gdot
real(rk)                :: u,alpha,fp,r0,v0s
real(rk)                :: a,asq,en
real(rk)                :: dm,ec,es,esq,xkep
real(rk)                :: fchk,s,c

!...  Executable code 

!...  Set dt = dt0 to be sure timestep is not altered while solving
!...  for new coords.
   dt = dt0
   iflg = 0_ik
   r0 = sqrt(dot_product(x0(1:3),x0(1:3)))
   v0s = dot_product(vx0(1:3),vx0(1:3))
   u = dot_product(x0(1:3),vx0(1:3))
   alpha = 2.0_rk*mu/r0-v0s

   if (alpha.gt.0.0_rk) then
      a = mu/alpha
      asq = a**2
      en = sqrt(mu/(a*asq))
      ec = 1.0-r0/a
      es = u/(en*asq)
      esq = ec**2+es**2
      dm = dt*en-int(dt*en/TWOPI)*TWOPI
      dt = dm/en
      if ( (dm**2.le.0.16_rk) .and. (esq.le.0.36_rk) ) then
         if (esq*dm**2 .lt. 0.0016_rk) then
            call drift_kepmd(dm,es,ec,xkep,s,c)
            fchk = (xkep - ec*s +es*(1.0-c) - dm)

            if(fchk**2 .gt. DANBYB) then
               iflg = 1_ik
               return
            endif

            fp = 1.0 - ec*c + es*s
            f = (a/r0) * (c-1.0) + 1.0
            g = dt + (s-xkep)/en
            fdot = - (a/(r0*fp))*en*s
            gdot = (c-1.0)/fp + 1.0

            x(:) = x0(:)*f + vx0(:)*g
            vx(:) = x0(:)*fdot + vx0(:)*gdot

            x0(:) = x(:)
            vx0(:) = vx(:)
            iflg = 0_ik
            return
         endif
      endif
   endif
       
   call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

   if (iflg .eq.0) then
        f = 1.0 - (mu/r0)*c2
        g = dt - mu*c3
        fdot = -(mu/(fp*r0))*c1
        gdot = 1.0 - (mu/fp)*c2

        x(:) = x0(:)*f + vx0(:)*g
        vx(:) = x0(:)*fdot + vx0(:)*gdot

        x0(:) = x(:)
        vx0(:) = vx(:)
   endif

return
end subroutine drift_dan
