!*************************************************************************
!                            DISCARD_MASS_PERI.F
!*************************************************************************
! This subroutine checks to see if a partical should be discarded because
! of its perihelion distance gets too small
!             Input:
!                 time           ==>  current time (real scalar)
!                 nbod           ==>  number of test bodies (int scalar)
!                 iecnt          ==>  Number of encounters (int*2 array)
!                 mass           ==>  mass of bodies (real array)
!                 xh,yh,zh       ==>   part position in helio coord 
!                                      (real arrays)
!                 vxh,vyh,vzh    ==>   part vel in helio coord 
!                                      (real arrays)
!                 qmin           ==>  Smallest perihelion distance 
!                                      (real scalar)
!                 iwhy           ==>  status of the object
!                                      (integer array)
!             Output:
!                 iwhy           ==>  status of the object
!                                      (integer array)
!                 isperi         ==> = 0 if tp went through peri
!                                    =-1 if tp pre peri
!                                    = 1 if tp post peri
!                                         (integer array)
! Remarks: Based on discard_peri
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision: 

subroutine discard_mass_peri5p(time,nbod,iecnt,mass,xh,                &
                               vxh,qmin,iwhy,isperi)
use swift_mod
use util_interface
implicit none

integer(ik), intent(in)    :: nbod,iecnt(:)
real(rk), intent(in)       :: mass(:),time,qmin
real(rk), intent(in)       :: xh(:,:),vxh(:,:)

integer(ik), intent(inout) :: iwhy(:)
integer(ik), intent(inout) :: isperi(:)

integer(ik)                :: i
integer(ik),save           :: i1st = 0_ik
real(rk)                   :: peri(NTPMAX)
logical(ik)                :: lperi(NTPMAX)

!...  Executable code

   if(i1st.eq.0) then     ! if first time through, set things up
      call util_mass_peri(0_ik,nbod,xh,vxh,mass,isperi,peri,lperi)
      i1st = 1_ik
      return              !  <==== RETURN
   endif

   call util_mass_peri(1_ik,nbod,xh,vxh,mass,isperi,peri,lperi)

   do i=2,nbod
      if ( (isperi(i).eq.0) .and. (iecnt(i).eq.0) ) then
         if (peri(i).le.qmin) then
            write(*,*) 'Particle',i,' perihelion distance too',        &
                       ' small at t=',time
            iwhy(i) = -4_ik
         endif
      endif
   enddo

return
end subroutine discard_mass_peri5p
