!*************************************************************************
!                            SYMBA5P_GROUP.F
!*************************************************************************
!             Input:
!                 ielst         ==>  list of ecnounters (2D integer array)
!                 iecnt         ==>  Number of encounters (int*2 array)
!                 i_ie          ==>  i_ie as in ielst(1,ielc)
!                 j_ie          ==>  j_ie as in ielst(2,ielc)
!             Output:
!                 indlst        ==>  list of independent pairs (integer array)
!                 indc          ==>  Number of independent pairs (integer)
!                 grpie(i,j)    ==>  lists of pairs in groups (2D integer array)
!                                    stores the index ie in ielst for pair no. i 
!                                    in group no. j
!                 grppc(i)      ==>  Number of pairs in group i (integer array)
!                 grpc          ==>  Number of groups (integer)
!             Internal:
!                 tempc         ==>  Counter of no. of groups that i_ie and j_ie 
!                                    are already in (integer)
!                 temp(2)       ==>  list of group(s) that i_ie and j_ie already in (integer)
!                 ie_t          ==>  The temporary pair no. (integer)
! Remarks:
! Authors:  CH Lau
! Date:    06/08/20

subroutine symba5p_group(ielst,ielc,i_ie,j_ie,grpie,grppc,grpc)
implicit none
use swift_mod
use sybam5p_mod
use util_interface

integer(ik), intent(in)  :: ielst(:,:),ielc

integer, intent(inout)   :: grpie(:,:),grppc(:),grpc

integer(ik)              :: i,j,i_ie,j_ie,ie_t
integer(ik)              :: temp(2),tempc

!...  Executable code
!...  group the encounter pairs
   tempc = 0_ik
   temp = 0_ik
!... see if i_ie or j_ie is in any group already
   do j=1,grpc
      do i=1,grppc(j)
      ie_t = grpie(i,j)
      if( ((i_ie.eq.ielst(1,ie_t)).or.
      &    (i_ie.eq.ielst(2,ie_t))) .or.
      &   ((j_ie.eq.ielst(1,ie_t)).or.
      &    (j_ie.eq.ielst(2,ie_t))) )then
         if (tempc.eq.0) then
            tempc = 1_ik
            temp(1) = j
            exit
         else
            tempc = 2_ik
            temp(2) = j
            exit
         endif
      endif
      enddo
      if (tempc.eq.2) then
         exit
      endif
   enddo
!... if i_ie or j_ie is not in any group yet, create and put this pair in a new group
   if (tempc.eq.0) then
      grpc = grpc+1
      if (grpc.gt.GRPNMAX) then
         write(*,*) 'ERROR: Group list filled.'
         write(*,*) 'STOPPING'
         call util_exit(1)
      endif
      grppc(grpc) = 1_ik
      grpie(1,grpc) = ielc
!... if i_ie and/or j_ie is already in any exisitng group(s), put this pair in the lower index group
   else
      grppc(temp(1)) = grppc(temp(1)) + 1
      if (grppc(temp(1)).ge.GRPMAX) then
         write(*,*) 'ERROR: Group member list filled.'
         write(*,*) 'STOPPING. Pair count=',grppc(temp(1))
         call util_exit(1)
      endif
      grpie(grppc(temp(1)),temp(1)) = ielc
   endif
!... if i_ie and j_ie already are in two different groups, merge larger index to lower index group 
   if (tempc.eq.2) then
      do i=1,grppc(temp(2))
         grppc(temp(1)) = grppc(temp(1))+1
         if (grppc(temp(1)).ge.GRPMAX) then
            write(*,*) 'ERROR: Group member list filled.'
            write(*,*) 'STOPPING. Pair count=',grppc(temp(1))
            call util_exit(1)
      	endif
         grpie(grppc(temp(1)),temp(1)) = grpie(i,temp(2))
         grpie(i,temp(2)) = 0_ik
      enddo
      grppc(temp(2)) = 0_ik
!... move one group up
      do j=(temp(2)+1),grpc
         do i=1,grppc(j)
            grpie(i,j-1) = grpie(i,j)
         enddo
         grppc(j-1) = grppc(j)
      enddo
      grpc = grpc-1
   endif

return
end subroutine symba5p_group
