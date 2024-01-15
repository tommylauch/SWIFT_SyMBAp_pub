c*************************************************************************
c                          ANAL_ENERGY_WRITE_SYMBAP.F
c*************************************************************************
c Writes the energy of the total system (massive bodies) wrt time.
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh            ==>  current position in helio coord 
c                               (real array)
c            vxh           ==>  current velocity in helio coord 
c                               (real array)
c            iu            ==>  unit to write to (int scalar)
c            fopenstat     ==>  The status flag for the open 
c                                statements of the output files.  
c                                      (character*80)
c            eoff          ==> An energy offset that is added to the energy
c                                      (real*8 scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/4/93
c Last revision: 12/27/96

subroutine anal_energy_write(t,nbod,mass,j2rp2,j4rp4,xh,vxh,           &
                             iu,fopenstat,eoff)
implicit none
use swift_mod

integer(ik), intent(in)        :: nbod,iu
real(rk), intent(in)           :: mass(:),t,j2rp2,j4rp4,eoff
real(rk), intent(in)           :: xh(:,:),vxh(:,:)
character(len = :), intent(in) :: fopenstat

integer(ik), save              :: i1st = 0_ik                           ! = 0 first time through; =1  after
real(rk)                       :: energy,eltot(3),ke,pot

c...  Executable code 

c Compute and print initial ke,pot,energy and ang. mom.
call anal_energy(nbod,mass,j2rp2,j4rp4,xh,vxh,ke,pot,energy,eltot)

energy = energy+eoff

call io_energy_write(i1st,t,energy,eltot,iu,fopenstat)

if(i1st.eq.0) i1st = 1_ik

return
end subroutine anal_energy_write

