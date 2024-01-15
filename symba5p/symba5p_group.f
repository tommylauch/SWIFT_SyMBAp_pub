c*************************************************************************
c                            SYMBA5P_GROUP.F
c*************************************************************************
c
c             Input:
c                 ielst         ==>  list of ecnounters (2D integer array)
c                 iecnt         ==>  Number of encounters (int*2 array)
c                 i_ie          ==>  i_ie as in ielst(1,ielc)
c                 j_ie          ==>  j_ie as in ielst(2,ielc)
c             Output:
c                 indlst        ==>  list of independent pairs (integer array)
c                 indc          ==>  Number of independent pairs (integer)
c                 grpie(i,j)    ==>  lists of pairs in groups (2D integer array)
c                                    stores the index ie in ielst for pair no. i 
c                                    in group no. j
c                 grppc(i)      ==>  Number of pairs in group i (integer array)
c                 grpc          ==>  Number of groups (integer)
c             Internal:
c                 tempc         ==>  Counter of no. of groups that i_ie and j_ie 
c                                    are already in (integer)
c                 temp(2)       ==>  list of group(s) that i_ie and j_ie already in (integer)
c                 ie_t          ==>  The temporary pair no. (integer)
c Remarks:
c Authors:  CH Lau
c Date:    06/08/20

      subroutine symba5p_group(ielst,ielc,i_ie,j_ie,grpie,grppc,grpc)
      include '../swift.inc'
      include '../symba5/symba5.inc'
      include 'symba5p.inc'

c...  Inputs Only: 
      integer ielst(:,:),ielc,i_ie,j_ie
c...  Outputs only
      integer grpie(:,:),grppc(:),grpc
c...  Internals
      integer i,j,ie_t
      integer temp(2),tempc

c...  Executable code
c...  group the encounter pairs
          tempc = 0
          temp = 0
c... see if i_ie or j_ie is in any group already
          do j=1,grpc
             do i=1,grppc(j)
                ie_t = grpie(i,j)
                if(((i_ie.eq.ielst(1,ie_t)).or.
     &             (i_ie.eq.ielst(2,ie_t))).or.
     &             ((j_ie.eq.ielst(1,ie_t)).or.
     &             (j_ie.eq.ielst(2,ie_t))))then
                   if (tempc.eq.0) then
                       tempc = 1
                       temp(1) = j
                       exit
                   else
                       tempc = 2
                       temp(2) = j
                       exit
                   endif
                endif
             enddo
             if (tempc.eq.2) then
                exit
             endif
          enddo
c... if i_ie or j_ie is not in any group yet, create and put this pair in a new group
          if (tempc.eq.0) then
             grpc = grpc + 1
             if (grpc.gt.GRPNMAX) then
                write(*,*) 'ERROR: Group list filled.'
                write(*,*) 'STOPPING'
                call util_exit(1)
             endif
             grppc(grpc) = 1
             grpie(1,grpc) = ielc
c... if i_ie and/or j_ie is already in any exisitng group(s), put this pair in the lower index group
          else
             grppc(temp(1)) = grppc(temp(1)) + 1
             if (grppc(temp(1)).ge.GRPMAX) then
                write(*,*) 'ERROR: Group member list filled.'
                write(*,*) 'STOPPING. Pair count=',grppc(temp(1))
                call util_exit(1)
             endif
             grpie(grppc(temp(1)),temp(1)) = ielc
          endif
c... if i_ie and j_ie already are in two different groups, merge larger index to lower index group 
          if (tempc.eq.2) then
             do i=1,grppc(temp(2))
                grppc(temp(1)) = grppc(temp(1)) + 1
                if (grppc(temp(1)).ge.GRPMAX) then
                   write(*,*) 'ERROR: Group member list filled.'
                   write(*,*) 'STOPPING. Pair count=',grppc(temp(1))
                   call util_exit(1)
             	endif
                grpie(grppc(temp(1)),temp(1)) = grpie(i,temp(2))
                grpie(i,temp(2)) = 0
             enddo
             grppc(temp(2)) = 0
c... move one group up
             do j=(temp(2)+1),grpc
                do i=1,grppc(j)
                   grpie(i,j-1) = grpie(i,j)
                enddo
                grppc(j-1) = grppc(j)
             enddo
             grpc = grpc - 1
          endif
      return
      end ! symba5p_group.f
c-----------------------------------------------------------
