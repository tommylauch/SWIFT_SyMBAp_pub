c*************************************************************************
c                            IO_SPLITNAME
c*************************************************************************
c splits the directory from the filename in a string
c
c             Input:
c                 oname        ==> string with the full path (character*80)
c   
c             Output:
c                 dirname     ==> string with the path (character*80)
c                 ldir        ==> length of dirname (integer scalar)
c                 filename     ==> string with the file name (character*80)
c                 lfile        ==> length of filename (integer scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:   3/19/97 
c Last revision: 

subroutine io_splitname(oname,dirname,ldir,filename,lfile)
implicit none
use swift_mod
use io_mod

character(len = :), intent(in)  :: oname

integer(ik), intent(out)        ::ldir,lfile
character(len = :), intent(out) :: dirname,filename

integer(ik)                     :: i,il,is

c...  Executable code 

c... Find the last character
il = 0
do i=1,80
   if (oname(i:i).eq.' ') then
      il = i-1
      exit
   endif
enddo

if (il.eq.0_ik) il = 80

c... Find the last '/'
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

