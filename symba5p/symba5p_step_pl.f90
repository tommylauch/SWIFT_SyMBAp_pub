!*************************************************************************
!                           SYMBA5P_STEP_PL.F
!*************************************************************************
!            Input:
!                i1st          ==>  = 0 if first step; = 1 not (int scalar)
!                time          ==>  current time (real scalar)
!                nbod          ==>  number of massive bodies (int scalar)
!                nbodm         ==>  location of the last massie body 
!                                   (int scalar)
!                mass          ==>  mass of bodies (real array)
!                j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                    (real scalars)
!                xh            ==>  initial position in helio coord
!                                   (real arrays)
!                vxh           ==>  initial velocity in helio coord
!                                   (real arrays)
!                dt            ==>  time step
!                lclose        ==> .true. --> marge particles if they
!                                   get too close. Read in that 
!                                   distance in io_init_pl
!                                     (logical*2 scalar)
!                rpl           ==>  physical size of a planet.
!                                   (real array)
!                eoff          ==>  Energy offset (real scalar)
!                rhill         ==>  size of planet's hills sphere
!                                   (real array)
!                mtiny         ==>  Small mass  (real array)
!            Output:
!                xh            ==>  final position in helio coord
!                                      (real arrays)
!                vxh           ==>  final velocity in helio coord
!                                      (real arrays)
!                rpl           ==>  Recalculated physical size of a planet.
!                                   if merger happened (real array)
!                nbod          ==>  Recalculated number of massive bodies 
!                                   if merger happened (int scalar)
!                mass          ==>  Recalculated mass of bodies 
!                                   if merger happened (real array)
!                isenc         ==>  0 --> No encounter during last dt
!                                   1 --> There was encounters
!                                    (integer scalar)
!                mergelst      ==>  list of mergers (int array)
!                mergecnt      ==>  count of mergers (int array)
!                iecnt         ==>  Number of encounters (int*2 array)
!                eoff          ==>  Energy offset (real scalar)
!                rhill         ==>  size of planet's hills sphere
!                                   (real array)
!                grpie(i,j)    ==>  lists of pairs in groups (2D integer array)
!                                   stores the index ie in ielst for pair no. i 
!                                   in group no. j
!                grppc(i)      ==>  Number of pairs in group i (integer array)
!                grpc          ==>  Number of groups (integer)
! Remarks: Based on symba2_step_pl.f 
! Authors:  Hal Levison
! Date:    11/27/97
! Last revision: 

subroutine symba5p_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,j4rp4,      &
                           xh,vxh,dt,lclose,rpl,isenc,                 &
                           mergelst,mergecnt,iecnt,eoff,rhill,mtiny)
implicit none
use swift_mod
use sybam5p_mod
use util_interface
use symba5p_interface, except_this_one => symba5p_step_pl

integer(ik), intent(in)  :: nbod,i1st,nbodm
real(rk), intent(in)     :: mass(:),dt,time,j2rp2,j4rp4,mtiny
logical(ik), intent(in)  :: lclose

real(rk), intent(inout)  :: xh(:,:),vxh(:,:)
real(rk), intent(inout)  :: rpl(:),eoff,rhill(:)

integer(ik), intent(out) :: isenc
integer(ik), intent(out) :: iecnt(:),ielev(:)
integer(ik), intent(out) :: mergelst(2,:),mergecnt

integer(ik)              :: i,j,ieflg,irec
logical(ik)              :: svdotr            ! Not used in the routine
integer(ik)              :: ielst(2,NENMAX),ielc
integer(ik)              :: grpie(GRPMAX,GRPNMAX),grppc(GRPNMAX),grpc

!...  Executable code 

   iecnt = 0_ik
   ielev = -1_ik
   isenc = 0_ik
!...  check for encounters
   irec = 0_ik
   ielc = 0_ik
   grpie = 0_ik
   grpc = 0_ik
   grppc = 0_ik
!$OMP PARALLEL DEFAULT (NONE)
!$OMP& PRIVATE(i,j,ieflg,svdotr)
!$OMP& SHARED(rhill,nbod,nbodm,mass,xh,vxh,
!$OMP& dt,irec,iecnt,ielev,ielc,ielst,isenc,grpie,grppc,grpc)
!$OMP DO COLLAPSE(2)
   do j=2,nbodm
      do i=j+1,nbod
         ieflg = 0_ik
         call symba5p_chk(rhill,nbod,i,j,mass,xh,vxh,dt,irec,ieflg,svdotr)
         if (ieflg.ne.0) then
!$OMP CRITICAL (ENC)
            isenc = 1_ik
            iecnt(i) = iecnt(i)+1
            iecnt(j) = iecnt(j)+1
            ielev(i) = 0_ik
            ielev(j) = 0_ik
            ielc = ielc+1
            if (ielc.gt.NENMAX) then
               write(*,*) 'ERROR: encounter matrix is filled.'
               write(*,*) 'STOPPING'
               call util_exit(1)
            endif
            ielst(1,ielc) = i
            ielst(2,ielc) = j
            call symba5p_group(ielst,ielc,i,j,grpie,grppc,grpc)
!$OMP END CRITICAL (ENC)
         endif
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!...  do a step
   if (isenc.eq.0) then
      call symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,        &
                              xh,vxh,dt)
      mergecnt = 0_ik
   else
      call symba5p_step_interp(time,iecnt,ielev,nbod,nbodm,mass,       &
           rhill,j2rp2,j4rp4,lclose,rpl,xh,vxh,dt,mergelst,mergecnt,   &
           eoff,ielc,ielst,mtiny,grpie,grppc,grpc)
      i1st = 0_ik
   endif

return
end subroutine symba5p_step_pl
