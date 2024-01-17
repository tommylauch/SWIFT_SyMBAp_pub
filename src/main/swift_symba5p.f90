program swift_symba5p
use swift_mod
use anal_interface
use util_interface
use io_interface
use symba5p_interface
use discard_interface
implicit none

real(rk)           :: mass(NTPMAX),xh(3,NTPMAX),vxh(3,NTPMAX)
real(rk)           :: xht(3,1),vxht(3,1)             ! Dummy for the io
real(rk)           :: j2rp2,j4rp4
integer(ik)        :: nbod,i1st,nbodm,nbodo,ntp,istat(1)

integer(ik)        :: threads,th_low,th_max
integer(ik)        :: iflgchk,iub,iuj,iud,iue,ium

real(rk)           :: t0,tstop,dt,dtout,dtdump
real(rk)           :: t,tout,tdump,tfrac,eoff
real(rk)           :: rpl(NTPMAX),rhill(NTPMAX)

real(rk)           :: rmin,rmax,rmaxu,qmin,mtiny
real(rk)           :: ke,pot,energy,eltot(3)
logical(ik)        :: lclose 
integer(ik)        :: isenc,ihills
integer(ik)        :: mergelst(2,NTPMAX),mergecnt,iecnt(NTPMAX)

character(len=:), allocatable :: outfile,inparfile,inplfile,fopenstat

!...  Executable code
   ntp = 0_ik

!...  print version number
   call util_version

! Get data for the run and the test particles
   write(*,*) 'Enter name of parameter data file : '
   read(*,'(a)') inparfile
   call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,iflgchk,      &
                      rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

! clean up discard_mass.out if needed
   if ((fopenstat(1:6).eq.'append').or.(fopenstat(1:6).eq.'APPEND')) then
      call io_discard_cleanup(t0)
   endif

! Prompt and read name of planet data file
   write(*,*) ' '
   write(*,*) 'Enter name of planet data file : '
   read(*,'(a)') inplfile
   call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,vxh,rpl,rhill, &
                   j2rp2,j4rp4)

   write(*,*) 'Enter the smallest mass to self gravitate :'
   read(*,*) mtiny
   write(*,*) ' mtiny = ',mtiny

! Get threads usage parameter:
   write(*,*) 'Enter lower limit of particle-to-thread ratio : '
   read(*,'(i6)') th_low

   write(*,*) 'Max. no. of threads to be used : '
   read(*,'(i6)') th_max

   threads = max(min(nbod/th_low,th_max),1)
   call omp_set_num_threads(threads)
   write(*,*) 'No. of threads: ',threads

! Initialize initial time and times for first output and first dump
   t = t0
   tout = t0 + dtout
   tdump = t0 + dtdump

   iub = 20_ik
   iuj = 30_ik
   iud = 40_ik
   iue = 60_ik
   ium = 21_ik

!...    Do the initial io write
   call io_write_frame(t0,nbod,ntp,mass,xh,vxh,xht,vxht,istat,         &
                       outfile,iub,fopenstat)
   call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)
!...  must initize discard io routine
   if (btest(iflgchk,4)) then ! bit 4 is set
      call io_discard_mass(0,t,0,mass(1),rpl(1),xh(1:3,1),vxh(1:3,1),  &
                           iud,-1,fopenstat)
   endif

!...  Calculate the location of the last massive particle
   call symba5p_nbodm(nbod,mass,mtiny,nbodm)

!...  set up energy write stuff
   if (btest(iflgchk,2)) then ! bit 2 is set
      eoff = 0.0_rk
      call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,vxh,          &
                             iue,fopenstat,eoff)
      call anal_energy_discard5(1,nbod,nbodm,mass,j2rp2,j4rp4,         &
                                xh,vxh,ke,pot,energy,eltot)
   else
      call anal_energy_discard5(-1,nbod,nbodm,mass,j2rp2,              &
                                j4rp4,xh,vxh,ke,pot,energy,eltot)
   endif

   ihills = 0_ik
   i1st = 0_ik

   write(*,*) ' ************** MAIN LOOP ****************** '
   do while ( (t.le.tstop).and.(nbod.gt.1) )
      call symba5p_step_pl(i1st,t,nbod,nbodm,mass,j2rp2,j4rp4,         &
                           xh,vxh,dt,lclose,rpl,isenc,                 &
                           mergelst,mergecnt,iecnt,eoff,rhill,mtiny)
      t = t+dt

      if (btest(iflgchk,4)) then     ! bit 4 set
         nbodo = nbod
         call discard_massive5p(t,dt,nbod,mass,xh,vxh,rmin,rmax,       &
              rmaxu,qmin,lclose,rpl,rhill,isenc,mergelst,mergecnt,     &
              iecnt,eoff,i1st)
         if (nbodo.ne.nbod) then
            call symba5p_nbodm(nbod,mass,mtiny,nbodm)
            if ( (nbod/threads.lt.th_low).and.(threads.gt.1) ) then
                ! change no. of threads
               threads = max(nbod/th_low,1)
               call omp_set_num_threads(threads)
               write(*,*) 'No. of threads decreased to ',threads
            endif
         endif
      endif

      if (t.ge.tout) then            ! output orb. elements
         call io_write_frame(t,nbod,ntp,mass,xh,vxh,xht,vxht,istat,    &
                             outfile,iub,fopenstat)
         call io_write_mass(t,nbod,mass,outfile,ium,fopenstat)
      tout = tout+dtout
      endif

      if (t.ge.tdump) then           ! do a dump
         tfrac = (t-t0)/(tstop-t0)
         write(*,fmt_prog) t,tfrac,nbod
         call io_dump_pl('dump_pl.dat',nbod,mass,xh,vxh,               &
                         lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
         call io_dump_param('dump_param.dat',t,tstop,dt,dtout,         &
              dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
         tdump = tdump+dtdump
         if (btest(iflgchk,2))  then ! bit 2 set
            call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,xh,vxh,     &
                                   iue,fopenstat,eoff)
         endif
      endif
   enddo

! Do a final dump for possible resumption later 
   call io_dump_pl('dump_pl.dat',nbod,mass,xh,vxh,                     &
                   lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
   call io_dump_param('dump_param.dat',t,tstop,dt,dtout,               &
        dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
        
   call util_exit(0)

end program swift_symba5p
