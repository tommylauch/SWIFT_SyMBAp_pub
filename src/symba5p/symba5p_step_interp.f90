!*************************************************************************
!                             SYMBA5P_STEP_INTERP.F
!*************************************************************************
!              Input:
!                  time          ==> Current time (real scalar)
!                  iecnt         ==>  The number of objects that each planet 
!                                     is encountering (int*2 array)
!                  ielev         ==>  The level that this particle should go
!                                              (int*2 array)
!                  nbod          ==>  number of massive bodies (int scalar)
!                  nbodm         ==>  Location of last massive body(int scalar)
!                  mass          ==>  mass of bodies (real array)
!                  rhill         ==>  Radius of hill sphere (real array)
!                  j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                      (real scalars)
!                  xh,yh,zh      ==>  initial position in helio coord 
!                                     (real arrays)
!                  vxh,vyh,vzh   ==>  initial velocity in helio coord 
!                                     (real arrays)
!                  dt            ==>  time step
!                  lclose        ==> .true. --> marge particles if they
!                                     get too close. Read in that 
!                                     distance in io_init_pl
!                                       (logical*2 scalar)
!                  rpl           ==>  physical size of a planet.
!                                     (real array)
!                  eoff          ==>  Energy offset (real scalar)
!                 iel!            ==>  number of encounters (integer scalar)
!                 ielst          ==>  list of ecnounters (2D integer array)
!              Output:
!                  xh,yh,zh      ==>  final position in helio coord 
!                                        (real arrays)
!                  vxh,vyh,vzh   ==>  final velocity in helio coord 
!                                        (real arrays)
!                  rpl           ==>  Recalculated physical size of a planet.
!                                     if merger happened (real array)
!                  nbod          ==>  Recalculated number of massive bodies 
!                                     if merger happened (int scalar)
!                  nbodm         ==>  Location of last massive body(int scalar)
!                  mass          ==>  Recalculated mass of bodies 
!                                     if merger happened (real array)
!                  mergelst      ==>  list of mergers (int array)
!                  mergecnt      ==>  count of mergers (int array)
!                  eoff          ==>  Energy offset (real scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:    11/21/96
! Last revision: 5/13/99

subroutine symba5p_step_interp(time,iecnt,ielev,nbod,nbodm,mass,       &
      rhill,j2rp2,j4rp4,lclose,rpl,xh,vxh,dt,mergelst,mergecnt,        &
      eoff,ielc,ielst,grpie,grppc,grpc)
use swift_mod
use symba5p_mod
use coord_interface
use helio_interface
use mvs_interface
use symba5p_interface, except_this_one => symba5p_step_interp
implicit none

real(rk), intent(in)       :: dt,j2rp2,j4rp4,time
integer(ik), intent(in)    :: nbod,nbodm,iecnt(:)
integer(ik), intent(in)    :: grpie(:,:),grpc
logical(ik), intent(in)    :: lclose

integer(ik), intent(inout) :: ielst(:,:),ielc,ielev(:),grppc(:)
real(rk), intent(inout)    :: xh(:,:),vxh(:,:),rpl(:),eoff
real(rk), intent(inout)    :: mass(:),rhill(:)

integer(ik), intent(out)   :: mergelst(:,:),mergecnt

integer(ik)                :: irec,ilevl(NTPMAX),i,j,k,ip2
integer(ik)                :: gplst(2,GRPMAX)
real(rk)                   :: dth
real(rk), save             :: axh(3,NTPMAX),vxb(3,NTPMAX)
real(rk)                   :: ptxb(3),ptxe(3)                           ! Not used here
logical(ik)                :: svdotr(NENMAX)                            ! Used by symba_step_recur

!...  Executable code 
   dth = 0.5_rk*dt

!...  Convert vel to bery to jacobi coords
   call coord_vh2b(nbod,mass,vxh,vxb)

!...  Do the linear drift due to momentum of the Sun
   call helio_lindrift(nbod,mass,vxb,dth,xh,ptxb)

!...  Get the accelerations in helio frame. For each object
!...     only include those guys that it is not encountering with. 
   call symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,axh,ielc,ielst)

!...  Apply a heliocentric kick for a half dt 
   call kickvh(nbod,vxb,axh,dth)

!..   Do a recursion step for full dt for particles not in close encounter
   irec = -1_ik
   call symba5p_helio_drift(nbod,ielev,irec,mass,xh,vxb,dt)
   irec = 0_ik
   do i=2,nbod
      ilevl(i) = 0_ik
   enddo
   mergecnt = 0_ik
!------------------------------
!$OMP PARALLEL DEFAULT (NONE)                                          &
!$OMP PRIVATE(i,j,svdotr,gplst)                                        &
!$OMP SHARED(ielst,time,nbod,nbodm,mass,irec,eoff,ielev,rhill,xh,vxb,  &
!$OMP lclose,rpl,mergelst,mergecnt,dt,grpc,grppc,grpie,iecnt,ilevl)
!$OMP DO SCHEDULE(GUIDED)
!...  Guided scheduling as groups have different sizes and go to different levels
   do j=1,grpc
      do i=1,grppc(j) !cheat recur step with gplst and grppc(j)
          gplst(1,i) = ielst(1,grpie(i,j))
          gplst(2,i) = ielst(2,grpie(i,j))
      enddo
   call symba5p_step_recur(time,nbod,nbodm,mass,irec,ilevl,iecnt,      &
        ielev,rhill,xh,vxb,lclose,rpl,mergelst,mergecnt,dt,eoff,       &
        svdotr,grppc(j),gplst)
   enddo
!$OMP END DO
!$OMP END PARALLEL
!..   Remove any encounters with ip2s in the encounter matrix
   do i=1,mergecnt
       ip2 = mergelst(2,i)
       j = 1_ik
       do while(j.le.ielc)
         if ( (ielst(1,j).eq.ip2) .or. (ielst(2,j).eq.ip2) ) then
            do k=j+1,ielc
               ielst(1,k-1) = ielst(1,k)
               ielst(2,k-1) = ielst(2,k)
            enddo
            ielc = ielc-1
         else
            j = j+1
         endif
       enddo
   enddo
!------------------------------

!...  Get the accelerations in helio frame. For each object
!...     only include those guys that it is not encountering with. 
   call symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,axh,ielc,ielst)

!...  Apply a heliocentric kick for a half dt 
   call kickvh(nbod,vxb,axh,dth)

!...  Do the linear drift due to momentum of the Sun
   call helio_lindrift(nbod,mass,vxb,dth,xh,ptxe)

!...  convert back to helio velocities
   call coord_vb2h(nbod,mass,vxb,vxh)

return
end subroutine symba5p_step_interp
