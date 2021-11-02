c*************************************************************************
c                            IO_OPEN_FXDR.F
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c open files using the fxdr routines
c
c             Input:
c                 fname           ==>  file name (character*80)
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*1)
c                 lflg            ==>  if .TRUE., then FXDR routines return even if there is
c                                      an I/O error.  If .FALSE., then routines halt on 
c                                      I/O error (logical scalar)
c
c             Output:
c                 iu              ==>  unit number (integer scalar)
c                 ierr            ==>  output from iostat
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/3/99
c Last revision: 

      subroutine io_open_fxdr(fname,fopenstat,lflg,iu,ierr)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu
      logical lflg
      character*(*) fname
      character*1 fopenstat

c...  Outputs: 
      integer ierr

c...  Internals:
      integer initxdr

c----
c...  Executable code 

#ifdef  _FXDR_AVAIL

      iu = initxdr(fname, fopenstat, lflg)

      if(iu.gt.0) then
         ierr = 0
      else
         ierr = iu
      endif

#else

      write(*,*) ' SWIFT ERROR: in io_open_fxdr: '
      write(*,*) '     This machine does not support FXDR'
      call util_exit(1)

#endif

      return
      end   ! io_open
