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

      subroutine symba5p_tune(tune,threads,start)

      include '../swift.inc'
      include '../symba5/symba5.inc'

! input
      integer start
! input and output
      integer tune,threads
! internal
      integer*8 sc_start,sc_end,sc_elapse,sc_elapseo
      integer step,threadso
      parameter (step=5000)
      save sc_start,sc_elapseo,threadso
c-----
c...  Executable code 

      if ((start.eq.1).and.(tune.eq.1)) then
         sc_elapseo = 0
         tune = 2
      endif

      if ((start.eq.1).and.(tune.eq.2)) then
         threadso = threads
         call system_clock(sc_start)
      endif

      if (start.eq.0) then
         if (tune.eq.step+2) then
            tune = 2
            call system_clock(sc_end)
            sc_elapse = sc_end-sc_start
            if ((sc_elapseo.eq.0) .or. (sc_elapse.lt.sc_elapseo)) then
               sc_elapseo = sc_elapse
               threadso = threads
               threads = threads-1
               if (threads .eq. 1) then
                  write(*,*) 'Tune done'
                  tune = 0
                  sc_elapseo = 0
               endif
            else
               threads = threadso
               write(*,*) 'Tune done'
               tune = 0
               sc_elapseo = 0
            endif
            call omp_set_num_threads(threads)
            write(*,*) 'No. of threads set to ',threads
         else
            tune = tune + 1
         endif
      endif
      return
      end                       ! symba5p_tune
c------------------------------------------------------


