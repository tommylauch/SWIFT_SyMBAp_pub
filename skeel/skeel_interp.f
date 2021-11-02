c*************************************************************************
c                            SKEEL_INTERP.F
c*************************************************************************
c This subroutine interpolates between two kepler orbits.
c
c             Input:
c                 msun                 ==>  mass of sun (real sclar)
c                 xbeg,ybeg,zbeg      ==>  initial planet position in helio 
c                                            (real scalars)
c                 vxbeg,vybeg,vzbeg   ==>  initial planet vel in bery
c                                            (real scalars)
c                 xend,yend,zend      ==>  final planet position in helio 
c                                            (real scalars)
c                 vxend,vyend,vzend   ==>  final planet position in bery
c                                            (real scalars)
c                 dti                 ==>  small time step (real sclar)
c             Output:
c                 xpl,ypl,zpl         ==>  position of planet wrt time 
c                                          for inner region
c                                            (real arrays)
c                 vxpl,vypl,vzpl      ==>  velcity of planet wrt time 
c                                          for inner region
c                                            (real arrays)
c
c
c Remarks: Based on rmvs_interp.f
c Authors:  Hal Levison 
c Date:    9/24/96
c Last revision: 3/18/97

      subroutine skeel_interp(msun,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,xpl,ypl,zpl,vxpl,vypl,
     &     vzpl,dti)

      include '../swift.inc'
      include 'skeel.inc'

c...  Inputs Only: 
      real*8 dti,msun
      real*8 xbeg,ybeg,zbeg
      real*8 vxbeg,vybeg,vzbeg
      real*8 xend,yend,zend
      real*8 vxend,vyend,vzend

c...  Outputs:
      real*8 xpl(0:NTENC),ypl(0:NTENC),zpl(0:NTENC)
      real*8 vxpl(0:NTENC),vypl(0:NTENC),vzpl(0:NTENC)
      
c...  Internals
      real*8 dt
      real*8 xc1,yc1,zc1,vxc1,vyc1,vzc1
      integer ib,iflg

c----
c...  Executable code 

      dt = dti*float(NTENC)

c...  move the end positions to beginning
      xc1 = xbeg
      yc1 = ybeg
      zc1 = zbeg
      vxc1 = vxbeg
      vyc1 = vybeg
      vzc1 = vzbeg

      xpl(0) = xbeg
      ypl(0) = ybeg
      zpl(0) = zbeg
      vxpl(0) = vxbeg
      vypl(0) = vybeg
      vzpl(0) = vzbeg

      do ib = 1,NTENC-1
            
         call drift_one(msun,xc1,yc1,zc1,vxc1,vyc1,vzc1,dti,iflg)
         if(iflg.ne.0) then
            write(*,*) ' Planet is lost in skeel_interp(2) !!!'
            write(*,*) msun,dti,ib
            write(*,*) xc1,yc1,zc1
            write(*,*) vxc1,vyc1,vzc1
            write(*,*) ' STOPPING '
            call util_exit(1)
         endif

         xpl(ib) = xc1
         ypl(ib) = yc1
         zpl(ib) = zc1
         vxpl(ib) = vxc1
         vypl(ib) = vyc1
         vzpl(ib) = vzc1

      enddo

      xpl(NTENC) = xend
      ypl(NTENC) = yend
      zpl(NTENC) = zend
      vxpl(NTENC) = vxend
      vypl(NTENC) = vyend
      vzpl(NTENC) = vzend

      return
      end      ! skeel_interp.f
c-----------------------------------------------------------------------

