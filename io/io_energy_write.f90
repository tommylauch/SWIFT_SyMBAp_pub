c*************************************************************************
c                          IO_ENERGY_WRITE.F
c*************************************************************************
c Does the write for anal_jacobi_write
c
c      Input:
c            i1st           ==>  =0 if first write, =1 if not (int scalar)
c            t              ==>  current time (real scalar)
c            energy         ==>  Total energy
c            eltot          ==>  components of total angular momentum
c                               (real array)
c            iu             ==>  unit to write to
c            fopenstat      ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/21/94
c Last revision: 3/4/94

subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
implicit none
use swift_mod
use io_mod

integer(ik), intent(in)        :: iu,i1st
real(rk), intent(in)           :: t,energy,eltot(:)
character(len = :), intent(in) :: fopenstat

integer(ik)                    :: ierr

c...  Executable code 

if(i1st.eq.0) then
   call io_open(iu,'energy.out',fopenstat,'FORMATTED',ierr)
   if(ierr.ne.0) then
      write(*,*) ' SWIFT ERROR: in anal_energy_write '
      write(*,*) '     Could not open energy.out '
      call util_exit(1)
   endif
else
   call io_open(iu,'energy.out','append','FORMATTED',ierr)
endif

write(iu,'1x,1p1e12.5,4(2x,1p1e23.16)') t,energy,eltot
close(iu)

return
end subroutine io_energy_write
