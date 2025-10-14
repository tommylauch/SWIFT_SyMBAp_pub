c*************************************************************************
c                            symba5p_tune.F
c*************************************************************************
c This subroutine tune the number of threads according to time taken
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

      subroutine symba5p_tune(tune,threads,threads_max,start)

      include '../swift.inc'
      include '../symba5/symba5.inc'

! input
      integer start,threads_max
! input and output
      integer tune,threads
! internal
      integer*8 sc_start,sc_end,sc_elapse,sc_elapseo
      integer ierr,step,threadso,threads_file
      character*80 fn
      parameter (fn = 'threads.dat')
      parameter (step=10000)
      save sc_start,sc_elapseo,threadso
c-----
c...  Executable code 

      if (start.eq.1) then
         if (tune.eq.1) then
            open(unit=10, file=fn, status='old', iostat=ierr)
            if (ierr.eq.0) then
               read(10,*) threads_file
               threads = max(threads,threads_file)
               close(10)
            else
               threads = threads_max
               write(*,*) fn,' not found, start tuning from ',threads
            endif
            call omp_set_num_threads(threads)
            write(*,*) 'Tunning, no. of threads set to ',threads
            sc_elapseo = 0
            threadso = threads
            tune = 2
         else if (tune.eq.2) then
            call system_clock(sc_start)
         endif
      else
         if (tune.eq.step+2) then
            call system_clock(sc_end)
            tune = 2
            sc_elapse = sc_end-sc_start
            if ((sc_elapseo.eq.0) .or. (sc_elapse.lt.sc_elapseo)) then
               sc_elapseo = sc_elapse
               threadso = threads
               threads = threads-1
               if (threads .eq. 0) then
                  threads = 1
                  tune = 0
               endif
            else
               threads = threadso
               tune = 0
            endif
            call omp_set_num_threads(threads)
            write(*,*) 'Tunning, no. of threads set to ',threads
         else
            tune = tune + 1
         endif
         if (tune.eq.0) then
            call omp_set_num_threads(threads)
            write(*,*) 'Tune done, no. of threads set to ',threads
            open(unit=10, file=fn, status='replace')
            sc_elapseo = 0
            write(10,*) threads
            close(10)
         endif
      endif

      return
      end                       ! symba5p_tune
c------------------------------------------------------


