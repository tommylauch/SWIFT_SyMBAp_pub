!*************************************************************************
!                            IO_SPLITNAME
!*************************************************************************
! splits the directory from the filename in a string
!             Input:
!                 oname        ==> string with the full path (character*80)
!   
!             Output:
!                 dirname     ==> string with the path (character*80)
!                 ldir        ==> length of dirname (integer scalar)
!                 filename     ==> string with the file name (character*80)
!                 lfile        ==> length of filename (integer scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:   3/19/97 
! Last revision: 

subroutine io_splitname(oname,dirname,ldir,filename,lfile)
implicit none
use swift_mod

character(len = :), intent(in)  :: oname

integer(ik), intent(out)        :: ldir,lfile
character(len = :), intent(out) :: dirname,filename

integer(ik)                     :: i,il,is

!...  Executable code 

!... Find the last character
   il = 0_ik
   do i=1,80
      if (oname(i:i).eq.' ') then
         il = i-1
         exit
      endif
   enddo

   if (il.eq.0_ik) il = 80_ik

!... Find the last '/'
   is = 0_ik
   do i=1,il
      if (oname(i:i).eq.'/') is = i
   enddo

   if (is.eq.0_ik) then          ! there is no path
      do i=1,il
         filename(i:i) = oname(i:i)
      enddo
      lfile = il
      write(dirname,'./')
      ldir = 2_ik
   else
      do i=1,is
         dirname(i:i) = oname(i:i)
      enddo
      ldir = is
      do i=is+1,il
         filename(i-is:i-is) = oname(i:i)
      enddo
      lfile = il-is
   endif

return
end subroutine io_splitname

