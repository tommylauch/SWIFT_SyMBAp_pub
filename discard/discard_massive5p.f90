!*************************************************************************
!                            DISCARD_MASSIVE5P.F
!*************************************************************************
! This subroutine checks to see if a massive body should be discarded or
! merged.
!             Input:
!                 time          ==>  current time (real scalar)
!                 dt            ==>  time step  (real scalar)
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 xh,yh,zh      ==>   position in helio coord 
!                                    (real arrays)
!                 vxh,vyh,vzh   ==>   pl vel in helio coord 
!                                    (real arrays)
!                 rmin,rmax      ==>  maximum and min distance from Sun
!                                     if <0  then don't check
!                                        (real scalar)
!                 rmaxu          ==>  maximum distance from Sun in not bound
!                                     if <0  then don't check
!                                        (real scalar)
!                  qmin          ==> Smallest perihelion distance
!                                      if <0  then don't check
!                                          (real scalar)
!                 lclose        ==> .true. --> marge particles if they
!                                    get too close. Read in that 
!                                    distance in io_init_pl
!                                      (logical*2 scalar)
!                 rpl           ==>  physical size of a planet.
!                                    (real array)
!                 rhill         ==>  size of a planet's hill's sphere.
!                                    (real array)
!                 isenc         ==>  0 --> No encounter during last dt
!                                    1 --> There was encounters
!                                     (integer scalar)
!                 eoff          ==> Amount of energy lost due to discards
!                                          (real scalar)
!                 mergelst      ==>  list of mergers (int array)
!                 mergecnt      ==>  count of mergers (int array)
!                 iecnt         ==>  Number of encounters (int*2 array)
!                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
!             Output:
!                 nbod          ==>  recalculated number of massive bodies 
!                                       (int scalar)
!                 mass          ==>  recalculated mass of bodies (real array)
!                 xh,yh,zh      ==>  recalculated position in helio coord 
!                                    (real arrays)
!                 vxh,vyh,vzh   ==>  recalculated pl vel in helio coord 
!                                    (real arrays)
!                 rpl           ==> recalculated physical sizes of a planet.
!                                    (real array)
!                 rhill         ==>  reordered size of planet's hill's sphere.
!                                    (real array)
!                 eoff          ==> Updated amount of energy lost from discards
!                                          (real scalar)
!                 i1st          ==>  set to 0 if reordered (int scalar)
! Remarks:
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision: 5/13/99

subroutine discard_massive5p(time,dt,nbod,mass,xh,vxh,rmin,rmax,rmaxu, &
          qmin,lclose,rpl,rhill,isenc,mergelst,mergecnt,iecnt,eoff,i1st)
implicit none
use swift_mod
use discard_interface, except_this_one => discard_massive5p
use coord_interface
use anal_interface
use io_interface

real(rk), intent(in)       :: time,dt
integer(ik), intent(in)    :: isenc
real(rk), intent(in)       :: rmin,rmax,rmaxu,qmin
logical(ik), intent(in)    :: lclose
integer(ik), intent(in)    :: mergelst(:,:),mergecnt,iecnt(:)

integer(ik), intent(inout) :: nbod,i1st
real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:)
real(rk), intent(inout)    :: eoff,rpl(:),rhill(:)

integer(ik)                :: iwhy(NTPMAX),i,iu,iflg,i1,i2,j
integer(ik), save          :: isperih(NTPMAX)
real(rk)                   :: xb(3,NTPMAX),vxb(3,NTPMAX)
real(rk)                   :: rmin2,rmax2,rmaxu2,energy
real(rk)                   :: ei,ef,ke,pot,eltot(3),vdotr
real(rk)                   :: rh2,rb2,vb2,msys
logical(ik)                :: lrflg(NTPMAX)
character(len = :)         :: cdummy

!...  Executable code

!.... check for duplicate mergers
   do i=1,nbod
      lrflg(i) = .true.
   enddo
   do i=1,mergecnt
      i2 = mergelst(2,i)
      if (lrflg(i2)) then
         lrflg(i2) = .false.
      else
         mergelst(2,i) = -1_ik
      endif
   enddo

!.... take care of mergers
   do i=1,mergecnt
      i1 = mergelst(1,i)
      i2 = mergelst(2,i)
      vdotr = dot_product(xh(1:3,i1),vxh(1:3,i1))
      if (vdotr.gt.0.0_rk) then
         isperih(i1) = 1_ik
      else 
         isperih(i1) = -1_ik
      endif
      if (i2.gt.0) then
         call discard_mass_reorder5(i2,nbod,mass,xh,vxh,rpl,           &
                                    rhill,isperih)
         i1st = 0_ik
         do j=i+1,mergecnt
            if (mergelst(1,j).gt.i2) then
               mergelst(1,j) = mergelst(1,j)-1
            endif
            if (mergelst(2,j).gt.i2) then
               mergelst(2,j) = mergelst(2,j)-1
            endif
         enddo
      endif
   enddo

!...  set things up
   iwhy = 0_ik
!...  check for position
   if ( (rmin.ge.0.0) .or. (rmax.ge.0.0) .or. (rmaxu.ge.0.0) ) then
      rmin2 = rmin**2
      rmax2 = rmax**2
      rmaxu2 = rmaxu**2
      call coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)
      do i=2,nbod
         rh2 = dot_product(xh(1:3,i),xh(1:3,i))
         if ( (rmax.ge.0.0) .and. (rh2.gt.rmax2) ) then
            write(*,*) rmax2,rh2,i
            write(*,*) 'Particle',i,' too far from Sun at t=',time
            iwhy(i) = -3_ik
         endif
         if ( (rmin.ge.0.0) .and. (rh2.lt.rmin2) ) then
            write(*,*) 'Particle',i,' too close from Sun at t=',time
            iwhy(i) = 1_ik
         endif

         if ((iecnt(i).eq.0).and.(rmaxu.ge.0.0).and.(iwhy(i).eq.0)) then
            rb2 = dot_product(xb(1:3,i),xb(1:3,i))
            vb2 = dot_product(vxb(1:3,i),vxb(1:3,i))
            energy = 0.5_rk*vb2-msys/sqrt(rb2)
            if ( (energy.gt.0.0) .and. (rb2.gt.rmaxu2) ) then
               write(*,*) 'Particle',i,' is unbound and too far ',     &
                          'from barycenter at t=',time
               iwhy(i) = -2_ik
            endif
         endif
      enddo
   endif

!...  check perihelion distance
   if (qmin.ge.0.0) then
      call discard_mass_peri5p(time,nbod,iecnt,mass,xh,vxh,            &
                               vxh,qmin,iwhy,isperih)
   endif
   iu = 40_ik
   i = 2_ik
   iflg = 0_ik
   do while(i.le.nbod) 
      if (iwhy(i).ne.0) then
         if (iflg.eq.0) then
            iflg = 1_ik
            call anal_energy(nbod,mass,0.0_rk,0.0_rk,xh,               &
                             vxh,ke,pot,ei,eltot)
         endif
         call io_discard_mass(1_ik,time,i,mass(i),rpl(i),xh(1:3,i),    &
                              vxh(1:3,i),iu,iwhy(i),cdummy)
         do j=i,nbod-1
            iwhy(j) = iwhy(j+1)
         enddo
         i1st = 0_ik
         call discard_mass_reorder5(i,nbod,mass,xh,vxh,rpl,            &
                                    rhill,isperih)
      else
         i = i+1
      endif
   enddo
   if (iflg.ne.0) then
      call anal_energy(nbod,mass,0.0_rk,0.0_rk,xh,                     &
                       vxh,ke,pot,ef,eltot)
      eoff = ei-ef
   endif

return
end subroutine discard_massive5p
