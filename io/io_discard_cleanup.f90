!*************************************************************************
!                            IO_DISCARD_CLEANUP
!*************************************************************************
! clean up discard_mass.out upon resumption
! Authors:  Tommy CH Lau
! Date:    14/03/2022

subroutine io_discard_cleanup(t0)
implicit none
use swift_mod

real(rk), intent(in) :: t0,time

integer(ik)          :: i,ierr,iu,it,ic,iw,iwhy,flg
character(len = :)   :: line

!...  Executable code 

   iu = 40_ik
   it = 41_ik
   ic = 42_ik
   flg = 0_ik

   open(iu,file='discard_mass.out',status='old',form='formatted',iostat=ierr)
   if(ierr.ne.0) then  
      write(*,*) ' discard_mass.out not found, no cleanup needed'
      return
   endif
   open(it,file='discard_mass.temp',status='new',form='formatted')
   open(ic,file='discard_mass.cleanup',access='append',form='formatted')

   do
      read(iu,*,iostat=ierr) time,iwhy
      if (ierr.ne.0) exit
      if (time.lt.t0) then
         iw = it
      else
         iw = ic
         flg = 1
      endif
      write(iw,'1x,1p1e23.16,1x,i4') time,iwhy
      if(iwhy.eq.2) then
         do i=1,9
            read(iu,'(A)') line
            write(iw,'(A)') line
         enddo
      else
         do i=1,3
            read(iu,'(A)') line
            write(iw,'(A)') line
         enddo
      endif
   enddo

   if (flg.eq.0) then
      close(iu)
      close(it,status='delete')
      close(ic,status='delete')
      write(*,*) ' No cleanup needed for discard_mass.out'
   else
      close(iu,status='delete')
      call rename('discard_mass.temp','discard_mass.out')
      close(it)
      close(ic)
      write(*,*) ' Cleanup done for discard_mass.out'
   endif

return
end subroutine io_discard_cleanup

