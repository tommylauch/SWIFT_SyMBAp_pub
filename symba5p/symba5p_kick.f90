!*************************************************************************
!                             SYMBA5P_KICK.F
!*************************************************************************
! Do a symba5 kick
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 irec          ==>  recursion level  (integer scalar)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 ielev         ==>  The level that this particle should go
!                                             (int*2 array)
!                 rhill         ==>  Hill sphere of planet (real Scalar)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb           ==>  initial velocity in bari coord 
!                                    (real arrays)
!                dt             ==>  timestep  (real scalar)
!                sgn            ==>  add or subtract the force (real scalar)
!                ielc           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
!            Output:
!                 vxb           ==>  final velocity in bari coord 
!                                    (real arrays)
! Remarks: Uses Man Hoi's force
! Authors:  Hal Levison 
! Date:   3/20/97
! Last revision: 3/3/10

subroutine symba5p_kick(nbod,mass,irec,iecnt,ielev,                    &
                        rhill,xh,vxb,dt,sgn,ielc,ielst)
implicit none
use swift_mod
use sybam5p_mod

integer(ik), intent(in) :: nbod,irec
real(rk), intent(in)    :: mass(:),dt,rhill(:),sgn
integer(ik), intent(in) :: iecnt(:),ielev(:),ielst(:,:),ielc
real(rk), intent(in)    :: xh(:,:)

real(rk), intent(inout) :: vxb(:,:)

real(rk)                :: dx(3),fac,ris,r
real(rk)                :: ri,rr,r2,faci,facj,ir3,rim1
integer(ik)             :: i,j,irm1,irecl,ie

!...  Executable code

   irm1 = irec-1
   if(sgn.lt.0.0_rk) then
      irecl = irec-1
   else
      irecl = irec
   endif

!...  calculate the accelerations
   do ie=1,ielc
      i = ielst(1,ie)
      j = ielst(2,ie)
      if((ielev(i).ge.irm1) .and. (ielev(j).ge.irm1) ) then
         ri = (rhill(i)+rhill(j))**2*RHSCALE**2*(RSHELL**(2*irecl))
         rim1 = ri*RSHELL**2
         
         dx(:) = xh(:,j)-xh(:,i)
         r2 = dot_product(dx,dx)
         ir3 = 1.0_rk/(r2*sqrt(r2))

         if (r2.lt.rim1) then
            fac = 0.0_rk
         else if (r2.lt.ri) then
            ris = sqrt(ri)
            r = sqrt(r2)
            rr = (ris-r)/(ris*(1.0_rk-RSHELL))
            fac = ir3 * (1.0_rk - (rr**4)*(35.0_rk -                   &
                         rr*(84.0_rk - rr*(70.0_rk - 20.0_rk*rr))))
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
end subroutine symba5p_kick
