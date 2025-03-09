c*************************************************************************
c                            SYMBA5P_THREAD.F
c*************************************************************************
c This subroutine changes the number of threads according to nbod
c
c             Input:
c                 nbodm         ==>  number of massive bodies (int scalar)
c                 nbod          ==>  total number of bodies (int scalar)
c                 th_low        ==>  limit of particle-to-thread ratio  (int scalar)
c                 th_max      ==>  initial position in helio coord 
c             Output:
c                 threads       ==> no. of thread to use (int scalar)
c
c Remarks: 
c Authors:  TCH Lau
c Date:   3/20/97
c Last revision: 
c

      subroutine symba5p_thread(nbod,nbodm,threads,th_low,th_max)

      include '../swift.inc'
      include '../symba5/symba5.inc'

! input
      integer nbod,nbodm,th_low,th_max
! input and output
      integer threads
! internal
      integer threadso

c-----
c...  Executable code 

      threadso = threads
      threads = min((nbodm+(nbod-nbodm)/2)/th_low+1,th_max)
      if (nbod .gt. 20) threads = max(threads,2)
      if (threads .ne. threadso) then
         call omp_set_num_threads(threads)
         write(*,*) 'No. of threads set to ',threads
      endif

      return
      end                       ! symba5p_thread
c------------------------------------------------------


