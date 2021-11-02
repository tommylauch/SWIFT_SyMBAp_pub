c*************************************************************************
c                             SKEEL_STEP_RECUR.F
c*************************************************************************
c                      THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c 
c
c             Input:
c                 msun          ==>  mass of Sun (real Scalar)
c                 mpl           ==>  mass of planet (real Scalar)
c                 rhill2        ==>  Hill sphere of planet (real Scalar)
c              xbeg,ybeg,zbeg   ==>  Initial Heliocentric position of planet 
c                                     (real Scalars)
c           vxbeg,vybeg,vzbeg   ==>  Initial Berycentric velocity of planet 
c                                     (real Scalars)
c              xend,yend,zend   ==>  Final Heliocentric position of planet 
c                                     (real Scalars)
c           vxend,vyend,vzend   ==>  Final Berycentric velocity of planet 
c                                     (real Scalars)
c                 xtp,ytp,ztp   ==>  Heliocentric position of TP
c                                     (real Scalars)
c              vxtp,vytp,vztp   ==>  Berycentric velocity of TP
c                                     (real Scalars)
c              istat1,istat2    ==>  first 2 istat flags (integer scalar)
c                        ireci  ==>  Input recursion level  (integer scalar)
c                        irecl  ==>  largest recursion level used (integer scalar)
c                        dt0    ==>  Global timestep  (real scalar)
c
c            Output:
c                 xtp,ytp,ztp   ==>  Heliocentric position of TP
c                                     (real Scalars)
c              vxtp,vytp,vztp   ==>  Berycentric velocity of TP
c                                     (real Scalars)
c              istat1,istat2    ==>  first 2 istat flags (integer scalar)
c                                    changes only if danby does not converge
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer scalar)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                    peri=QFLG if not good orbit or outside
c                                    the hills sphere of planet.
c                                         (real scalar)
c
c
c Remarks: 
c Authors:  Hal Levison 
c Date:   9/23/96
c Last revision: 2/13/01

#ifdef  _RECUR_SUB
      recursive subroutine 
#else
      subroutine 
#endif
     &     skeel_step_recur(msun,mpl,rhill2,xbeg,ybeg,zbeg,
     &         vxbeg,vybeg,vzbeg,xend,yend,zend,
     &         vxend,vyend,vzend,xtp,ytp,ztp,vxtp,vytp,vztp,
     &         istat1,istat2,ireci,irecl,dt0,isperi,peri)

      include '../swift.inc'
      include 'skeel.inc'

c...  Inputs Only: 
      real*8 msun,mpl,dt0,rhill2
      real*8 xbeg,ybeg,zbeg
      real*8 vxbeg,vybeg,vzbeg
      real*8 xend,yend,zend
      real*8 vxend,vyend,vzend

c...  Inputs & Outputs Only: 
      real*8 xtp,ytp,ztp
      real*8 vxtp,vytp,vztp
      real*8 peri
      integer ireci,istat1,istat2,irecl
      integer isperi

c...  Internals: 
      integer it,itm1,irecp
      real*8 dtl,dth,rim1,sgn,ri
      integer icflg,iflg
      real*8 xpl(0:NTENC),ypl(0:NTENC),zpl(0:NTENC)
      real*8 vxpl(0:NTENC),vypl(0:NTENC),vzpl(0:NTENC)
      real*8 xrel,yrel,zrel,vxrel,vyrel,vzrel,tperi,rrel2

c----
c...  Executable code 

      ri = rhill2*(RHSCALE*(RSHELL**ireci))**2

      dtl = dt0/(float(NTENC)**(ireci))
      dth = dtl/2.0

      irecl = max(irecl,ireci)

      if(ireci.eq.0) then

         irecp = ireci + 1

         call skeel_chk(mpl,xbeg,ybeg,zbeg,
     &        vxbeg,vybeg,vzbeg,
     &        xtp,ytp,ztp,vxtp,vytp,vztp,dtl,ri,icflg)

         sgn = 1.0d0
         call skeel_kick(mpl,xbeg,ybeg,zbeg,
     &        xtp,ytp,ztp,vxtp,vytp,vztp,ri,dth,sgn)

         if(icflg.eq.0) then
            call drift_one(msun,xtp,ytp,ztp,vxtp,vytp,vztp,dtl,iflg)
            if(iflg.ne.0) then
               istat1 = 1
               istat2 = -1
            endif
            xrel = xtp - xend
            yrel = ytp - yend
            zrel = ztp - zend
            rrel2 = xrel**2 + yrel**2 + zrel**2
            vxrel = vxtp - vxend
            vyrel = vytp - vyend
            vzrel = vztp - vzend
            call util_peri1(1,xrel,yrel,zrel,vxrel,vyrel,
     &          vzrel,mpl,isperi,peri,tperi)
            if(isperi.eq.0) then
               if( (abs(tperi).gt.(FACQDT*dtl)) .or. 
     &              (rrel2.gt.rhill2) ) then
c                  peri=QFLG
                  peri=sqrt(rrel2)
               endif
            endif
         else
            call skeel_step_recur(msun,mpl,rhill2,xbeg,
     &           ybeg,zbeg,vxbeg,vybeg,
     &           vzbeg,xend,yend,zend,vxend,
     &           vyend,vzend,xtp,ytp,ztp,vxtp,vytp,vztp,
     &           istat1,istat2,irecp,irecl,dt0,isperi,peri)
         endif

         call skeel_kick(mpl,xend,yend,zend,xtp,ytp,
     &        ztp,vxtp,vytp,vztp,ri,dth,sgn)

      else

         call skeel_interp(msun,xbeg,ybeg,zbeg,
     &        vxbeg,vybeg,vzbeg,xend,yend,zend,
     &        vxend,vyend,vzend,xpl,ypl,zpl,vxpl,vypl,vzpl,dtl)

         rim1 = rhill2*(RHSCALE*(RSHELL**(ireci-1)))**2
         
         irecp = ireci + 1
         do it=1,NTENC

            itm1 = it - 1
            call skeel_chk(mpl,xpl(itm1),ypl(itm1),zpl(itm1),
     &           vxpl(itm1),vypl(itm1),vzpl(itm1),
     &           xtp,ytp,ztp,vxtp,vytp,vztp,dtl,ri,icflg)

            sgn = 1.0d0
            call skeel_kick(mpl,xpl(itm1),ypl(itm1),zpl(itm1),xtp,
     &           ytp,ztp,vxtp,vytp,vztp,ri,dth,sgn)
            sgn = -1.0d0
            call skeel_kick(mpl,xpl(itm1),ypl(itm1),zpl(itm1),xtp,
     &           ytp,ztp,vxtp,vytp,vztp,rim1,dth,sgn)
            
            if(icflg.eq.0) then
               call drift_one(msun,xtp,ytp,ztp,vxtp,vytp,vztp,dtl,iflg)
               if(iflg.ne.0) then
                  istat1 = 1
                  istat2 = -1
               endif
               if(isperi.eq.-1) then
                  xrel = xtp - xpl(it)
                  yrel = ytp - ypl(it)
                  zrel = ztp - zpl(it)
                  rrel2 = xrel**2 + yrel**2 + zrel**2
                  vxrel = vxtp - vxpl(it)
                  vyrel = vytp - vypl(it)
                  vzrel = vztp - vzpl(it)
                  call util_peri1(1,xrel,yrel,zrel,vxrel,vyrel,
     &                 vzrel,mpl,isperi,peri,tperi)
                  if(isperi.eq.0) then
                     if( (abs(tperi).gt.(FACQDT*dtl)) .or. 
     &                    (rrel2.gt.rhill2) ) then
c                        peri=QFLG
                        peri=sqrt(rrel2)
                     endif
                  endif
               endif
            else
               call skeel_step_recur(msun,mpl,rhill2,xpl(itm1),
     &           ypl(itm1),zpl(itm1),vxpl(itm1),vypl(itm1),
     &           vzpl(itm1),xpl(it),ypl(it),zpl(it),vxpl(it),
     &           vypl(it),vzpl(it),xtp,ytp,ztp,vxtp,vytp,vztp,
     &              istat1,istat2,irecp,irecl,dt0,isperi,peri)
            endif

            sgn = 1.0d0
            call skeel_kick(mpl,xpl(it),ypl(it),zpl(it),xtp,
     &           ytp,ztp,vxtp,vytp,vztp,ri,dth,sgn)
            sgn = -1.0d0
            call skeel_kick(mpl,xpl(it),ypl(it),zpl(it),xtp,
     &           ytp,ztp,vxtp,vytp,vztp,rim1,dth,sgn)

         enddo

      endif

      return
      end     ! skeel_step_recur.f
c--------------------------------------------------------------
