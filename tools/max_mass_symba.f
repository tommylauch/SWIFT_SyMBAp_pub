c converts binary file to ascii file

      include 'swift.inc'

      real*8 mass(NTPMAX)
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      integer tot_inject
      integer iflgchk,iu,i,id,nleft,ium
      
      real*8 rpl(NTPMAX),rhill(NTPMAX)
      logical*2 lclose
      real*8 capm,j2rp2,j4rp4

      character*80 inplfile
      
      lclose = .true.
      iflgchk = 0
c Prompt and read name of planet data file
      read(*,999) inplfile
999   format(a)
      call io_init_pl_symba(inplfile,lclose,iflgchk,tot_inject,mass,xh,
     &       yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)
     
      mass(1) = 0.d0
      write(*,*) maxval(mass)/1.1922482116515945e-4
      
      return
      end
