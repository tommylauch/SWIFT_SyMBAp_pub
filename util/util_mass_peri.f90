!*************************************************************************
!                            UTIL_MASS_PERI.F
!*************************************************************************
! This subroutine determines whether peri of a planet has taken place
!             Input:
!                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
!                 nbod           ==>  number of bodies (int scalar)
!                 x,y,z          ==>  heliocentric position of planets
!                                       (real arrays)
!                 vx,vy,vz       ==>  heliocentric velcocities of planets
!                                       (real arrays)
!                 mass           ==>  mass of the bodies (real array)
!             Output:
!                 isperi         ==> = 0 if tp went through peri
!                                    =-1 if tp pre peri
!                                    = 1 if tp post peri
!                                         (integer array)
!                 peri           ==> set to pericenter dist. if isperi=0
!                                         (real array)
!                lperi           ==> set to .true. if isperi=0
!                                         (logical*2 array)
! Remarks: Based on util_peri.f
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision: 

subroutine util_mass_peri(iflg,nbod,x,vx,mass,isperi,peri,lperi)
implicit none
use swift_mod
use orbel_interface

integer(ik), intent(in)  :: nbod,iflg
real(rk), intent(in)     :: mass(:),x(:,:),vx(:,:),gm

real(rk), intent(out)    :: peri(:)
integer(ik), intent(out) :: isperi(:)
logical(ik), intent(out) :: lperi(:)

integer(ik)              :: i,ialpha
real(rk)                 :: vdotr,a,e

!...  Executable code 

   if (iflg.eq.0) then    ! are we just setting thing up?
      do i=2,nbod
         vdotr = dot_product(x(1:3,i),vx(1:3,i))
         if (vdotr .gt. 0.0_rk) then
            isperi(i) = 1_ik
         else 
            isperi(i) = -1_ik
         endif
      enddo
   else
      do i=2,nbod
         vdotr = dot_product(x(1:3,i),vx(1:3,i))
         if (isperi(i).eq.-1) then         ! was coming in
            if (vdotr .lt. 0.0_rk) then    ! still coming in
               isperi(i) = -1_ik
            else                           ! turned around
               isperi(i) = 0_ik
               lperi(i) = .true.
               gm = mass(1) + mass(i)
               call orbel_xv2aeq(x(1:3,i),vx(1:3,i),gm,ialpha,a,e,     &
                                 peri(i))
            endif
         else
            if (vdotr.lt.0.0_rk) then      ! coming in
               isperi(i) = -1_ik
            else
               isperi(i) = 1_ik            ! going out
            endif
         endif
      enddo
   endif

return
end subroutine util_mass_peri
