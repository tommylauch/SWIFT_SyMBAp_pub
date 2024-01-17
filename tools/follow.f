c converts binary file to ascii file

       include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri,sla_dranrm
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT),nskp,ic
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep,initxdr
	integer iflgchk,iu,nleft,i,id
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 apo,obar

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

        iu = 20

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an FXDR binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading an real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

        write(*,*) ' Input the particle number to follow '
        read(*,*) ifol
        write(*,*) ' Following particle ',ifol

	write(*,*) 'Input the print frequency'
	read(*,*) nskp

	if(btest(iflgchk,0))  then ! bit 0 is set
           call io_open_fxdr(outfile, 'r', .true., iu, ierr)
	else
           call io_open(iu,outfile,'old','UNFORMATTED',ierr)
	endif
        open(unit=7,file='follow.out')

        write(*,*) '1 2 3  4    5     6    7    8    9   10'
        write(*,*) 't,a,e,inc,capom,omega,capm,peri,apo,obar'

        tmax = t0
	ic = 0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif

             if(ierr.ne.0) then
		write(*,*) ' Stop while reading header ierr=',ierr
		goto 2
	     endif

             istep = 0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
		if(ierr.ne.0) then
		   write(*,*) ' Stop while reading planets ierr=',ierr
		   goto 2
		endif
                if(id.eq.ifol) then
                   istep = 1
                   obar = capom + omega
                   obar = obar/DEGRAD
		   obar = obar*dr
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
		   ic = ic + 1
		   if(mod(ic,nskp).eq.0) then
		      write(7,1000) t,a,e,inc,capom,omega,
     &                     capm,peri,apo,obar
 1000		      format(1x,e15.7,1x,f10.4,1x,f7.5,4(1x,f9.4),
     &                     2(1x,f10.4),1x,f9.4)
		   endif
                   tmax = t
                endif
             enddo

             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
		if(ierr.ne.0) then
		   write(*,*) ' Stop while reading planets ierr=',ierr
		   goto 2
		endif
                if(id.eq.ifol) then
                   istep = 1
                   obar = capom + omega
                   obar =  obar/DEGRAD
		   obar = obar*dr
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
		   ic = ic + 1
		   if(mod(ic,nskp).eq.0) then
		      write(7,1000) t,a,e,inc,capom,omega,
     &                     capm,peri,apo,obar
		   endif
                   tmax = t
                endif
             enddo
             if(istep.eq.0) goto 2     ! did not find particle this times step

        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

        stop
        end
