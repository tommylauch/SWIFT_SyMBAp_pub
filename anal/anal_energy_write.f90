!*************************************************************************
!                          ANAL_ENERGY_WRITE_SYMBAP.F
!*************************************************************************
! Writes the energy of the total system (massive bodies) wrt time.
!      Input:
!            t             ==>  current time
!            nbod          ==>  number of massive bodies (int scalar)
!            mass          ==>  mass of bodies (real array)
!            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
!            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
!            xh            ==>  current position in helio coord 
!                               (real array)
!            vxh           ==>  current velocity in helio coord 
!                               (real array)
!            iu            ==>  unit to write to (int scalar)
!            fopenstat     ==>  The status flag for the open 
!                                statements of the output files.  
!                                      (character*80)
!            eoff          ==> An energy offset that is added to the energy
!                                      (real*8 scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:    3/4/93
! Last revision: 12/27/96

subroutine anal_energy_write(t,nbod,mass,j2rp2,j4rp4,xh,vxh,           &
                             iu,fopenstat,eoff)
implicit none
use swift_mod
use anal_interface
use io_interface

integer(ik), intent(in)        :: nbod,iu
real(rk), intent(in)           :: mass(:),t,j2rp2,j4rp4,eoff
real(rk), intent(in)           :: xh(:,:),vxh(:,:)
character(len = :), intent(in) :: fopenstat

integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after
real(rk)                       :: energy,eltot(3),ke,pot

!...  Executable code 

! Compute and print initial ke,pot,energy and ang. mom.
   call anal_energy(nbod,mass,j2rp2,j4rp4,xh,vxh,ke,pot,energy,eltot)
   energy = energy+eoff
   call io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
   if (i1st.eq.0) i1st = 1_ik

return
end subroutine anal_energy_write

