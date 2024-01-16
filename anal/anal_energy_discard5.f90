!*************************************************************************
!                          ANAL_ENERGY_DISCARD5_SYMBAP.F
!*************************************************************************
! Calculates the energy of the total system (massive bodies) wrt time.
! returns the total energy of n objects by direct pairwise summation
! G = 1., and we allow diff. masses.  Also returns square of total ang. mom.
!      Input:
!            iflg          ==>  use to turn on/off energy calculation
!            t             ==>  current time
!            nbod          ==>  number of massive bodies (int scalar)
!            nbodm         ==>  Location of last massive body(int scalar)
!            mass          ==>  mass of bodies (real array)
!            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
!            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
!            xh            ==>  current position in heliocentric coord 
!                               (real array)
!            vxh           ==>  current velocity in heliocentric coord 
!                               (real array)
!      Output:
!            ke            ==>  kinetic energy
!            pot           ==>  potential energy
!            energy        ==>  Total energy
!            eltot         ==>  components of total angular momentum
!                               (real array)
! Remarks: Based on anal_energy
! Authors:  Hal Levison
! Date:  12/16/06
! Last revision:  

subroutine anal_energy_discard5(iflg,nbod,nbodm,mass,j2rp2,j4rp4,      &
                                xh,vxh,ke,pot,energy,eltot)
implicit none
use swift_mod
use anal_interface

integer(ik), intent(in) :: iflg,nbod,nbodm
real(rk), intent(in)    :: mass(:),j2rp2,j4rp4
real(rk), intent(in)    :: xh(:,:),vxh(:,:)

real(rk), intent(out)   :: energy,eltot(:),ke,pot

logical(ik), save       :: leuse = .true.

!...  Executable code
   if (iflg.lt.0) then
      leuse = .false.
      return
   else if (iflg.gt.0) then
      leuse = .true.
      return
   endif

   if (leuse) then
      call anal_energy_mtiny(nbod,nbodm,mass,j2rp2,j4rp4,xh,vxh,       &
                             ke,pot,energy,eltot)
   else
      ke = 0.0_rk
      pot = 0.0_rk
      energy = 0.0_rk
      eltot = 0.0_rk
   endif

return
end subroutine anal_energy_discard5
