!*************************************************************************
!                            HELIO_LINDRIFT.F
!*************************************************************************
! This subroutine takes a linear drift due to mometum of Sun
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 vxb           ==>  velocity in bary coord 
!                                    (real arrays)
!                 dt            ==>  time step
!                 xh            ==>  initial position in helio coord 
!                                       (real arrays)
!             Output:
!                 xh            ==>  final position in helio coord 
!                                       (real arrays)
!                 ptx           ==> momentum of sun: tp's need this   
!                                       (real scalars)
! Remarks: Bases on Martin's code h2.f
! Authors:  Hal Levison 
! Date:    11/14/96
! Last revision: 1/8/97

subroutine helio_lindrift(nbod,mass,vxb,dt,xh,ptx)
use swift_mod
implicit none

integer(ik), intent(in)  :: nbod
real(rk), intent(in)     :: mass(:),dt,vxb(:,:)

real(rk), intent(inout)  :: xh(:,:)

real(rk), intent(out)    :: ptx(:)

integer(ik)              :: n

!...  Executable code 

   ptx(:) = mass(2)*vxb(:,2)

   do n=3,nbod
      ptx(:) = ptx(:)+mass(n)*vxb(:,n)
   enddo

   ptx(:) = ptx(:)/mass(1)

   do n=2,nbod
      if(mass(n).ne.0.0_rk) then
         xh(:,n) = xh(:,n)+ptx(:)*dt
      endif
   enddo

return
end subroutine helio_lindrift
