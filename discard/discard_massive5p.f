c*************************************************************************
c                            DISCARD_MASSIVE5P.F
c*************************************************************************
c This subroutine checks to see if a massive body should be discarded or
c merged.
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rmin,rmax      ==>  maximum and min distance from Sun
c                                     if <0  then don't check
c                                        (real scalar)
c                 rmaxu          ==>  maximum distance from Sun in not bound
c                                     if <0  then don't check
c                                        (real scalar)
c                  qmin          ==> Smallest perihelion distance
c                                      if <0  then don't check
c                                          (real scalar)
c                 lclose        ==> .true. --> marge particles if they
c                                    get too close. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 rhill         ==>  size of a planet's hill's sphere.
c                                    (real array)
c                 isenc         ==>  0 --> No encounter during last dt
c                                    1 --> There was encounters
c                                     (integer scalar)
c                 eoff          ==> Amount of energy lost due to discards
c                                          (real scalar)
c                 mergelst      ==>  list of mergers (int array)
c                 mergecnt      ==>  count of mergers (int array)
c                 iecnt         ==>  Number of encounters (int*2 array)
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c             Output:
c                 nbod          ==>  recalculated number of massive bodies 
c                                       (int scalar)
c                 mass          ==>  recalculated mass of bodies (real array)
c                 xh,yh,zh      ==>  recalculated position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  recalculated pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==> recalculated physical sizes of a planet.
c                                    (real array)
c                 rhill         ==>  reordered size of planet's hill's sphere.
c                                    (real array)
c                 eoff          ==> Updated amount of energy lost from discards
c                                          (real scalar)
c                 i1st          ==>  set to 0 if reordered (int scalar)
c
c
c Remarks: 
c
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 5/13/99

      subroutine discard_massive5p(time,dt,nbod,mass,xh,vxh,rmin,rmax,
     &     rmaxu,qmin,lclose,rpl,rhill,isenc,mergelst,mergecnt,iecnt,
     &     eoff,i1st)

      include '../swift.inc'

c...  Inputs: 
      real*8 time,dt
      integer isenc
      real*8 rmin,rmax,rmaxu,qmin
      logical*2 lclose
      integer mergelst(2,NTPMAX),mergecnt
      integer iecnt(NTPMAX)

c...  Input and Output
      integer nbod,i1st
      real*8 mass(nbod),xh(3,nbod),vxh(3,nbod)
      real*8 eoff,rpl(nbod),rhill(nbod)

c...  internal
      integer iwhy(NTPMAX),i,iu,iflg,i1,i2,j,isperih(NTPMAX)
      real*8 xb(3,NTPMAX),vxb(3,NTPMAX)
      real*8 rmin2,rmax2,rmaxu2,energy
      real*8 ei,ef,ke,pot,eltot(3),vdotr
      real*8 rh2,rb2,vb2,msys
      logical*1 lrflg(NTPMAX)
      character*1 cdummy

      save isperih

c-----
c...  Executable code 


c.... check for duplicate mergers
      do i=1,nbod
         lrflg(i) = .true.
      enddo
      do i=1,mergecnt
         i2 = mergelst(2,i)
         if(lrflg(i2)) then
            lrflg(i2) = .false.
         else
            mergelst(2,i) = -1
         endif
      enddo

c.... take care of mergers
      do i=1,mergecnt
         i1 = mergelst(1,i)
         i2 = mergelst(2,i)
         vdotr =xh(1,i1)*vxh(1,i1)+xh(2,i1)*vxh(2,i1)+xh(3,i1)*vxh(3,i1)
         if (vdotr .gt. 0.d0) then
            isperih(i1) = 1
         else 
            isperih(i1) =-1
         endif
         if(i2.gt.0) then
            call discard_mass_reorder5_symbap(i2,nbod,mass,xh,vxh,rpl,
     &                                        rhill,isperih)
            i1st = 0
            do j=i+1,mergecnt
               if(mergelst(1,j).gt.i2) then
                  mergelst(1,j) = mergelst(1,j) - 1
               endif
               if(mergelst(2,j).gt.i2) then
                  mergelst(2,j) = mergelst(2,j) - 1
               endif
            enddo
         endif
      enddo

c...  set things up
      do i=1,nbod
         iwhy(i) = 0
      enddo

c...  check for position
      if( (rmin.ge.0.0) .or. (rmax.ge.0.0) .or. (rmaxu.ge.0.0) ) then
         rmin2 = rmin*rmin
         rmax2 = rmax*rmax
         rmaxu2 = rmaxu*rmaxu
         call coord_h2b_symbap(nbod,mass,xh,vxh,xb,vxb,msys)

         do i=2,nbod
            rh2 = xh(1,i)**2 + xh(2,i)**2 + xh(3,i)**2
            if( (rmax.ge.0.0) .and. (rh2.gt.rmax2) ) then
               write(*,*) rmax2,rh2,i
               write(*,*) 'Particle',i,' too far from Sun at t=',
     &              time
               iwhy(i) = -3
            endif
            if( (rmin.ge.0.0) .and. (rh2.lt.rmin2) ) then
               write(*,*) 'Particle',i,' too close from Sun at t=',
     &              time
               iwhy(i) = 1
            endif

            if((iecnt(i).eq.0).and.(rmaxu.ge.0.0).and.
     &           (iwhy(i).eq.0)) then
               rb2 = xb(1,i)**2 + xb(2,i)**2 + xb(3,i)**2
               vb2 = vxb(1,i)**2 + vxb(2,i)**2 + vxb(3,i)**2
               energy = 0.5*vb2 - msys/sqrt(rb2)
               if( (energy.gt.0.0) .and. (rb2.gt.rmaxu2) ) then
                  write(*,*) 'Particle',i,' is unbound and too far ',
     &                 'from barycenter at t=',time
                  iwhy(i) = -2
               endif
            endif
         enddo
      endif

c...  check perihelion distance
      if(qmin.ge.0.0) then
         call discard_mass_peri5p_symbap(time,nbod,iecnt,mass,xh,vxh,
     &                                   vxh,qmin,iwhy,isperih)
      endif

      iu = 40
      i = 2
      iflg = 0
      do while(i.le.nbod) 
         if(iwhy(i).ne.0) then
            if(iflg.eq.0) then
               iflg = 1
               call anal_energy_symbap(nbod,mass,0.0d0,0.0d0,xh,
     &           vxh,ke,pot,ei,eltot)
            endif
            call io_discard_mass_symbap(1,time,i,mass(i),rpl(i),xh(:,i),
     &           vxh(:,i),iu,iwhy(i),cdummy)
            do j=i,nbod-1
               iwhy(j) = iwhy(j+1)
            enddo
            i1st = 0
            call discard_mass_reorder5_symbap(i,nbod,mass,xh,vxh,rpl,
     &                                        rhill,isperih)
         else
            i = i + 1
         endif
      enddo


      if(iflg.ne.0) then
         call anal_energy_symbap(nbod,mass,0.0d0,0.0d0,xh,
     &        vxh,ke,pot,ef,eltot)
         eoff = ei - ef
      endif

      return
      end         ! discard_massive5p.f
c-----------------------------------------------------

