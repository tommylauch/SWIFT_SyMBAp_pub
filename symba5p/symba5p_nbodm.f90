!*************************************************************************
!                            SYMBA5_NBODM.F
!*************************************************************************
! Returns the location of the last massive body in the list
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 mtiny         ==>  Small mass  (real array)
!             Output:
!                 nbodm         ==>  location of the last massive body 
!                                    (int scalar)
! Remarks:  If all the objects are massive, then nbodm=nbod-1 so that
!           the do loops will have the correct limits.
! Authors:  Hal Levison
! Date:    3/20/97
! Last revision: 1/29/06

subroutine symba5p_nbodm(nbod,mass,mtiny,nbodm)
implicit none
use swift_mod
use sybam5p_mod

integer(ik), intent(in)  :: nbod
real(rk), intent(in)     :: mass(nbod),mtiny

integer(ik), intent(out) :: nbodm

integer(ik)              :: i

!...  Executable code

   if (mass(nbod).gt.mtiny) then
      nbodm = nbod-1
      write(*,*) '    Of ',nbod,' objects, all are massive. '
   else
      nbodm = 1_ik
      do i=1,nbod-1
         if (mass(i).gt.mtiny) nbodm = i
      enddo
      write(*,*) '    Of ',nbod,' objects, ',nbodm,' are massive. '
   endif

return
end subroutine symba5_nbodm
