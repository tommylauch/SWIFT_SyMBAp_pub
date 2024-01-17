! converts binary file to ascii file
program follow_symba5p
use swift_mod
use io_mod
implicit none

interface
   subroutine follow_plist(iflg,ifol,ifoln,plist,nbod,tg,tstop)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: iflg,ifol,nbod
   real(rk), intent(in)           :: tstop
   integer(ik), intent(inout)     :: plist(:)
   integer(ik), intent(out)       :: ifoln
   real(rk), intent(out)          :: tg
   end subroutine follow_plist
end interface

real(rk)               :: mass(NTPMAX),dr
real(rk)               :: xh(3,NTPMAX),vxh(3,NTPMAX)

integer(ik)            :: nbodm,nbod,nbod0,ierr,ifol,istep,tot_inject
integer(ik)            :: iflgchk,iu,i,id,nleft,ium
integer(ik)            :: io_read_hdr,io_read_line,io_read_mass

real(rk)               :: t0,tstop,dt,dtout,dtdump
real(rk)               :: t,tmax

real(rk)               :: rmin,rmax,rmaxu,qmin,rpl(NTPMAX),rhill(NTPMAX)
logical(ik)            :: lclose
real(rk)               :: a,e,inc,capom,omega,capm,j2rp2,j4rp4
real(rk)               :: peri,apo,tg

integer(ik)            :: plist(NTPMAX),ifoln

character(len=50)      :: outfile,inparfile,inplfile,fopenstat

! Get data for the run and the test particles
   write(*,*) 'Enter name of parameter data file : '
   read(*,'(a)') inparfile
   call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,iflgchk,      &
                      rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

! Prompt and read name of planet data file
   write(*,*) ' '
   write(*,*) 'Enter name of planet data file : '
   read(*,'(a)') inplfile
   call io_init_pl(inplfile,lclose,iflgchk,tot_inject,mass,xh,vxh,     &
                   rpl,rhill,j2rp2,j4rp4)

   iu = 20_ik
   ium = 30_ik
   dr = 180.0_rk/PI

   write(*,*) ' Reading an real*8 binary file '

   write(*,*) ' Input the particle number to follow '
   read(*,*) ifol
   ifol = abs(ifol)
   write(*,*) ' Following particle ',ifol

   open(unit=iu, file=outfile, status='old',form='unformatted')
   open(unit=ium, file='mass.'//outfile, status='old',form='unformatted')
   open(unit=7,file='follow_symba.out')

   nbod0 = tot_inject
   do i=1,nbod0
      plist(i) = i
   enddo
   call follow_plist(0_ik,ifol,ifoln,plist,nbod0,tg,tstop)

   ifoln = ifol

   write(*,*) '1  2 3 4  5    6     7    8    9   10  11 '
   write(*,*) 't,id,a,e,inc,capom,omega,capm,peri,apo, M '

   tmax = t0

   do
      ierr = io_read_hdr(iu,t,nbod,nleft) 
      ierr = io_read_mass(t,nbodm,mass,ium)
      if (ierr.ne.0_ik) exit
      if(nbodm.ne.nbod) then
         write(*,*) ' Error 1:',nbod,nbodm
         stop
      endif

      do while(t.ge.tg)
         call follow_plist(1_ik,ifol,ifoln,plist,nbod0,tg,tstop)
      enddo

      istep = 0_ik
      do i=2,nbod
         ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
         if (ierr.ne.0_ik) exit
         if (abs(id).eq.ifoln) then
            istep = 1_ik
            inc = inc*dr
            capom = capom*dr
            omega = omega*dr
            capm = capm*dr
            peri = a*(1.0_rk-e)
            apo = a*(1.0_rk+e)
            
            write(7,'(1x,e15.7,1x,i3,1x,f10.4,1x,f7.5,4(1x,f9.4),2(1x,f10.4),1e13.5)') &
            t,ifoln,a,e,inc,capom,omega,capm,peri,apo,mass(abs(id))/mass(1)
            
            tmax = t
         endif
      enddo
      if(ierr.ne.0_ik) exit
   enddo

        write(*,*) ' Tmax = ',tmax

stop
end program follow_symba5p
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine follow_plist(iflg,ifol,ifoln,plist,nbod,tg,tstop)
use swift_mod
implicit none
interface
   subroutine left_reorder(ig,im,nbod,plist)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: ig,im,nbod
   integer(ik), intent(inout)     :: plist(:)
   end subroutine left_reorder
end interface

integer(ik), intent(in)        :: iflg,ifol,nbod
real(rk), intent(in)           :: tstop

integer(ik), intent(inout)     :: plist(:)

integer(ik), intent(out)       :: ifoln
real(rk), intent(out)          :: tg

integer(ik), save              :: iwhy
integer(ik)                    :: ig,im,idum,i,ierr

   if (iflg.eq.0_ik) then
      open(2,file='discard_mass.out',status='old',iostat=ierr)
      if (ierr.ne.0_ik) then
         write(*,*) 'Could not open discard_mass.out'
         tg = 5.0_rk*tstop
         return            ! <====== NOTE 
      endif
      read(2,*,iostat=ierr) tg,iwhy
      if (ierr.ne.0_ik) then
         write(*,*) 'Could not read discard_mass.out'
         tg = 5.0_rk*tstop
         return            ! <====== NOTE 
      endif
      ifoln = ifol
      return               ! <====== NOTE 
   endif

   if (iwhy.eq.2) then
      read(2,*) idum,im
      read(2,fmt='(1x)')
      read(2,fmt='(1x)')
      read(2,*) idum,ig
      call left_reorder(ig,im,nbod,plist)
      do i=1,5
         read(2,fmt='(1x)')
      enddo
   else
      read(2,*) idum,ig
      im = -1
      call left_reorder(ig,im,nbod,plist)
      read(2,fmt='(1x)')
      read(2,fmt='(1x)')
   endif

   read(2,*,iostat=ierr) tg,iwhy
   if (ierr.ne.0) tg = 5.0_rk*tstop

   ifoln = plist(ifol)

return
end subroutine follow_plist
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine left_reorder(ig,im,nbod,plist)
use swift_mod
implicit none

integer(ik), intent(in)        :: ig,im,nbod
integer(ik), intent(inout)     :: plist(:)
integer(ik)                    :: i

   do i=1,nbod
      if(plist(i).eq.ig) then
         if(im.gt.0) then
            plist(i) = im
         else
            plist(i) = -1_ik
         endif
      endif
   enddo

   do i=1,nbod
      if(plist(i).gt.ig) plist(i) = plist(i)-1_ik
   enddo

return
end subroutine left_reorder
