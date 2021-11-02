c*************************************************************************
c                            IO_WRITE_MASS
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c NEW VERSION OF THIS, USES FXDR
c write out masses
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Hal Levison 
c Date:    11/2/99
c Last revision: 

      subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 mass(nbod),time
      character*80 oname,fopenstat

c...  Internals
      real*4 mass4(NTPMAX)
      real*4 ttmp
      integer ierr,i,ldir,lfile,initxdr,ixdrclose,ixdrreal,ixdrint
      integer ixdrrmat
      character*80 dirname,filename

      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st,dirname,filename,ldir,lfile

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

c...  if first time through open file
      if(i1st.eq.0) then
         call io_splitname(oname,dirname,ldir,filename,lfile)
         if( (fopenstat(1:6).eq.'append') .or. 
     &        (fopenstat(1:6).eq.'APPEND') ) then
            call io_open_fxdr(dirname(1:ldir)//'mass.'//
     &           filename(1:lfile),'a', .true., iu, ierr)
         else
            if( (fopenstat(1:3).eq.'new') .or. 
     &           (fopenstat(1:3).eq.'NEW') ) then
               call io_open_fxdr(dirname(1:ldir)//'mass.'//
     &              filename(1:lfile),'w', .true., iu, ierr)
               if(iu.ge.0) then
                  write(*,*) ' SWIFT ERROR: in io_write_mass: '
                  write(*,*) '     binary output file exists:'
                  call util_exit(1)
               endif
            endif
            call io_open_fxdr(dirname(1:ldir)//'mass.'//
     &           filename(1:lfile),'w', .true., iu, ierr)
         endif

         if(ierr.lt.0) then
           write(*,*) ' SWIFT ERROR: in io_write_mass: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
         call io_open_fxdr(dirname(1:ldir)//'mass.'//filename(1:lfile),
     &        'a', .true., iu, ierr)
         if(ierr.lt.0) then
           write(*,*) ' SWIFT ERROR: in io_write_mass: '
           write(*,*) '     Could not open binary output file'
           write(*,*) '      with append:'
           call util_exit(1)
         endif
      endif

      do i=1,nbod
         mass4(i) = mass(i)
      enddo
      ttmp = time

      ierr = ixdrreal( iu, ttmp )
      if(ierr.lt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_mass: '
         write(*,*) '     Could not write time'
         call util_exit(1)
      endif
      ierr = ixdrint( iu, nbod )
      if(ierr.lt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_mass: '
         write(*,*) '     Could not write nbod'
         call util_exit(1)
      endif
      ierr = ixdrrmat ( iu, nbod, mass4 )
      if(ierr.lt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_mass: '
         write(*,*) '     Could not write masses'
         call util_exit(1)
      endif

      ierr = ixdrclose(iu)
      if(ierr.lt.0) then
         write(*,*) ' SWIFT ERROR: in io_write_frame: '
         write(*,*) '     Could not close mass output file:'
         call util_exit(1)
      endif

#else

      write(*,*) ' SWIFT ERROR: in io_write_frame: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end      ! io_write_mass
c----------------------------------------------------------------------
