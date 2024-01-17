!*************************************************************************
!                            SYMBA5P_STEP_RECUR.F
!*************************************************************************
!                      THIS FILE MUST BE PRECOMPILED
!*************************************************************************
!             Input:
!                 t             ==>  time (real Scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 nbodm         ==>  Location of last massive body(int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 ireci         ==>  Input recursion level  (integer scalar)
!                 ilevl         ==>  largest recursion level used 
!                                    (integer array)
!                 iecnt         ==>  The number of objects that each planet 
!                                    is encountering (int*2 array)
!                 ielev         ==>  The level that this particle should go
!                                             (int*2 array)
!                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
!                                     (real scalars)
!                 rhill         ==>  Hill sphere of planet (real array)
!                 xh            ==>  initial position in helio coord 
!                                    (real arrays)
!                 vxb           ==>  initial velocity in bari coord 
!                                    (real arrays)
!                dt0            ==>  Global timestep  (real scalar)
!                lclose         ==> .true. --> marge particles if they
!                                    get too close. Read in that 
!                                    distance in io_init_pl
!                                      (logical*2 scalar)
!                rpl            ==>  physical size of a planet.
!                                    (real array)
!                eoff           ==>  Energy offset (real scalar)
!                iel!           ==>  number of encounters (integer scalar)
!                ielst          ==>  list of ecnounters (2D integer array)
!             Output:
!                 xh            ==>  final position in helio coord 
!                                       (real arrays)
!                 vxb           ==>  final velocity in bari coord 
!                                       (real arrays)
!             mergelst          ==>  list of mergers (int array)
!             mergecnt          ==>  count of mergers (int array)
!                 rpl           ==>  Recalculated physical size of a planet.
!                                    if merger happened (real array)
!                 mass          ==>  Recalculated mass of bodies 
!                                    if merger happened (real array)
!                eoff           ==>  Energy offset (real scalar)
!                svdotr         ==> vdotr relative flag
!                                   = .true. if i,j are receding
!                                   = .false is approaching
!                                     (2D logical*1 array)
!                                   Used internally, but only need 1 copy.
! Remarks: Based on SYMBA5_STEP_RECUR
! Authors:
! Date:

recursive subroutine symba5p_step_recur(t,nbod,nbodm,mass,ireci,ilevl, &
                     iecnt,ielev,rhill,xh,vxb,lclose,rpl,              &
                     mergelst,mergecnt,dt0,eoff,svdotr,ielc,ielst)
use swift_mod
use symba5p_mod
use symba5p_interface, except_this_one => symba5p_step_recur
implicit none

integer(ik), intent(in)    :: nbod,ireci,nbodm
real(rk), intent(in)       :: dt0,t
integer(ik), intent(in)    :: iecnt(:)
logical(ik), intent(in)    :: lclose

integer(ik), intent(inout) :: ilevl(:),ielst(:,:),ielc,ielev(:)
real(rk), intent(inout)    :: xh(:,:),vxb(:,:),eoff,rpl(:)
real(rk), intent(inout)    :: mass(:),rhill(:)
integer(ik), intent(inout) :: mergelst(:,:),mergecnt
logical(ik), intent(inout) :: svdotr(:)

integer(ik)                :: i,j,ie
integer(ik)                :: icflg,it,irecp,ieflg
real(rk)                   :: dtl,dth,sgn

!...  Executable code
   dtl = dt0/(float(NTENC)**(ireci))
   dth = dtl/2.0d0
   if ( (dtl/dt0) .le. TINY ) then
      write(*,*) ' Warning in SYMBA_STEP_RECUR: '
      write(*,*) '         Local timestep too small '
      write(*,*) '         Roundoff will be important!!!! '
   endif

   if (ireci.eq.0) then  ! Do we need to go deeper?
      irecp = ireci+1
      icflg = 0_ik
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)
         if ((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
            call symba5p_chk(rhill,i,j,xh,vxb,dtl,irecp,ieflg,svdotr(ie))
            if (ieflg.ne.0) then
               icflg = 1_ik
               ielev(i) = irecp
               ielev(j) = irecp
               ilevl(i) = max(irecp,ilevl(i))
               ilevl(j) = max(irecp,ilevl(j))
            endif
         endif
      enddo
      sgn = 1.0_rk
      call symba5p_kick(mass,irecp,iecnt,ielev,                        &
                        rhill,xh,vxb,dth,sgn,ielc,ielst)
      call symba5p_helio_drift_g(ielev,ireci,mass,xh,vxb,              &
                                 dtl,ielc,ielst)
      if (icflg.ne.0) then
         call symba5p_step_recur(t,nbod,nbodm,mass,irecp,ilevl,        &
              iecnt,ielev,rhill,xh,vxb,lclose,rpl,mergelst,mergecnt,   &
              dt0,eoff,svdotr,ielc,ielst)
      endif
      sgn = 1.0_rk
      call symba5p_kick(mass,irecp,iecnt,ielev,                        &
                        rhill,xh,vxb,dth,sgn,ielc,ielst)
      if ( lclose ) then       ! look for mergers
         do ie=1,ielc
            i = ielst(1,ie)
            j = ielst(2,ie)
            if ((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
               call symba5p_merge(t,dtl,i,j,mass,xh,vxb,svdotr(ie),    &
                    rpl,mergelst,mergecnt,rhill,eoff,ielc,ielst)
            endif
         enddo
      endif
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)
         if (ielev(i).eq.irecp) ielev(i) = ireci
         if (ielev(j).eq.irecp) ielev(j) = ireci
      enddo
   else
      irecp = ireci+1
      do it=1,NTENC                                                     ! Do we need to go deeper?
         icflg = 0_ik
         do ie=1,ielc
            i = ielst(1,ie)
            j = ielst(2,ie)
            if ((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
               call symba5p_chk(rhill,i,j,xh,vxb,dtl,irecp,ieflg,svdotr(ie))
               if (ieflg.ne.0) then
                  icflg = 1_ik
                  ielev(i) = irecp
                  ielev(j) = irecp
                  ilevl(i) = max(irecp,ilevl(i))
                  ilevl(j) = max(irecp,ilevl(j))
               endif
            endif
         enddo
         sgn = 1.0_rk
         call symba5p_kick(mass,irecp,iecnt,ielev,                     &
                           rhill,xh,vxb,dth,sgn,ielc,ielst)
         sgn = -1.0_rk
         call symba5p_kick(mass,irecp,iecnt,ielev,                     &
                           rhill,xh,vxb,dth,sgn,ielc,ielst)
         call symba5p_helio_drift_g(ielev,ireci,mass,xh,vxb,           &
                                    dtl,ielc,ielst)
         if (icflg.ne.0) then
            call symba5p_step_recur(t,nbod,nbodm,mass,irecp,ilevl,     &
                iecnt,ielev,rhill,xh,vxb,lclose,rpl,mergelst,          &
                mergecnt,dt0,eoff,svdotr,ielc,ielst)
         endif
         sgn = 1.0_rk
         call symba5p_kick(mass,irecp,iecnt,ielev,                     &
                           rhill,xh,vxb,dth,sgn,ielc,ielst)
         sgn = -1.0_rk
         call symba5p_kick(mass,irecp,iecnt,ielev,                     &
                           rhill,xh,vxb,dth,sgn,ielc,ielst)
         if (lclose) then ! look for mergers
            do ie=1,ielc
               i = ielst(1,ie)
               j = ielst(2,ie)
               if ((ielev(i).ge.ireci).and.(ielev(j).ge.ireci)) then
                  call symba5p_merge(t,dtl,i,j,mass,xh,vxb,svdotr(ie), &
                       rpl,mergelst,mergecnt,rhill,eoff,ielc,ielst)
               endif
            enddo
         endif
         do ie=1,ielc
            i = ielst(1,ie)
            j = ielst(2,ie)
            if (ielev(i).eq.irecp) ielev(i) = ireci
            if (ielev(j).eq.irecp) ielev(j) = ireci
         enddo
      enddo
   endif

return
end subroutine symba5p_step_recur
