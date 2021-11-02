c*************************************************************************
c                            IO_WRITE_FRAME
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c 
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  current position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  current velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  current part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  current velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c             Output
c                 iu              ==> unit number to write to.  
c                                     FXDR changes this
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      subroutine io_write_frame(time,nbod,ntp,mass,xh,yh,zh,vxh,
     &           vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,oname,
     &           iu,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,iu
      real*8 mass(nbod),time
      integer istat(ntp)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      character*(*) oname,fopenstat

c...  Internals
      integer i,id,ixdrclose
      integer ialpha, ierr
      real*8 a,e,inc,capom,omega,capm
      real*8 gm
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

c...  if first time through open file
      if(i1st.eq.0) then
         if( (fopenstat(1:6).eq.'append') .or. 
     &        (fopenstat(1:6).eq.'APPEND') ) then
            call io_open_fxdr(oname, 'a', .true., iu, ierr)
         else
            if( (fopenstat(1:3).eq.'new') .or. 
     &           (fopenstat(1:3).eq.'NEW') ) then
               call io_open_fxdr(oname, 'r', .true., iu, ierr)
               if(iu.ge.0) then
                  write(*,*) ' SWIFT ERROR: in io_write_frame: '
                  write(*,*) '     binary output file exists:'
                  call util_exit(1)
               endif
            endif
            call io_open_fxdr(oname, 'w', .true., iu, ierr)
         endif
         if(ierr.lt.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
         call io_open_fxdr(oname, 'a', .true., iu, ierr)
         if(ierr.lt.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file'
           write(*,*) '      with append:'
           call util_exit(1)
         endif
      endif

      call io_write_hdr(iu,time,nbod,ntp,istat)
      
c...  write out planets
      do i=2,nbod
         gm = mass(1)+mass(i)
         id = -1*i
 	 call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm,
     &          ialpha,a,e,inc,capom,omega,capm)
         call io_write_line(iu,id,a,e,inc,capom,omega,capm)
      enddo

c...  write out test particles
      gm = mass(1)
      do i=1,ntp
         if(istat(i).eq.0) then
            call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &          vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
            call io_write_line(iu,i,a,e,inc,capom,omega,capm)
         endif
      enddo

      ierr = ixdrclose(iu)
      if(ierr.lt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_frame: '
         write(*,*) '     Could not close binary output file:'
         call util_exit(1)
      endif

#else

      write(*,*) ' SWIFT ERROR: in io_write_frame: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end      ! io_write_frame
c----------------------------------------------------------------------
