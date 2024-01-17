c   Gererate initial position and velocity for test particles
c   All of them have the same e and i.  The user supplies a range in a. 
c   The tp are equally spaced in a.  The rest of the angles are chosen 
c   at random.

      include 'swift.inc'

      real*8 SMASSYR
      parameter(SMASSYR=TWOPI*TWOPI)

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 gm,a,e,inc,capom,omega,capm,amin,amax,da
      real*8 ran,rstat(NTPMAX,NSTATR)
      real*8 cimax,imax

      integer istat(NTPMAX,NSTAT)
      integer ntp,ialpha,i,j,iseed,iuflg

      real*4 ran2

      character*80 intpfile

 1    write(*,*) ' Units Menu: '
      write(*,*) '       0 ==> Solar masses, and AU '
      write(*,*) '       1 ==> AU, and Years '
      write(*,*) '       2 ==> other '
      read(*,*) iuflg
      if( (iuflg.ne.0) .and. (iuflg.ne.1)  .and. (iuflg.ne.2) ) goto 1

      if(iuflg.eq.0) then
         gm = 1.0d0
      else if(iuflg.eq.1) then
         gm = SMASSYR
      else 
         write(*,*) 'Input mass of the Sun'
         read(*,*) gm
      endif

      write(*,*) ' Input number particles '
      read(*,*) ntp
      write(*,*) ' Input ecc of particles '
      read(*,*) e
      write(*,*) ' Input max inc of particles (in deg) '
      read(*,*) imax
      write(*,*) ' Input amin and amax '
      read(*,*) amin,amax
      da = (amax-amin)/float(ntp-1)

c      write(*,*) ' Input iseed: large and odd '
c      read(*,*) iseed
      open(2,file='/home/hal/iseed.dat',status='old')
      read(2,*) iseed
      close(2)

      iseed = mod(abs(iseed),500000)
      iseed = -1 * (2*iseed + 1 )
      write(*,*) 'iseed = ',iseed

      a=amin
      cimax = cos(imax/DEGRAD)

      open(12,file='gen_ai.out')

      write(*,*) ' Test Paticles: '
      write(*,*) '   a      e      i    omega  capom    M '
      write(12,*) '   a      e      i    omega  capom    M '

      do i=1,ntp
         ialpha = -1
         inc = acos( 1.0d0 - (1.0d0-cimax)*ran2(iseed))
         capm = TWOPI*ran2(iseed)
         omega= TWOPI*ran2(iseed)
         capom = TWOPI*ran2(iseed)
         write(*,1000) a,e,inc,omega,capom,capm
         write(12,1000) a,e,inc,omega,capom,capm
 1000    format(6(1x,f6.3))

         call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &        xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i))
         do j=1,NSTAT
            istat(i,j) = 0
         enddo
         do j=1,NSTATR
            rstat(i,j) = 0.0d0
         enddo
         a = a + da
      enddo

      write(*,*) 'Enter name of test particle data file : '
      read(*,1001) intpfile
 1001 format(a)
      call io_dump_tp(intpfile,ntp,xht,yht,zht,
     &     vxht,vyht,vzht,istat,rstat)

      open(2,file='/home/hal/iseed.dat',status='old')
      write(2,*) iseed
      close(2)
      
      stop
      end    !  gen_a

c-------------------------------------------------------------------------
