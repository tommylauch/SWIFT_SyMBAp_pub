c************************************************************************
c                         IO_DUMP_PL_SYMBA.F
c************************************************************************
c Dumps the data for the Sun and planets 

c
c             Input:
c                 dplfile       ==>  Name of file to write to (character*80)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh            ==>  initial position in Helio coord 
c                                    (real arrays)
c                 vxh           ==>  initial position in Helio coord 
c                                    (real arrays)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
c                 rpl           ==>  physical size of planet
c                                    (real array)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c
c
c Remarks: Based on io_dump_pl.f
c Authors:  Hal Levison
c Date:    1/8/97
c Last revision: 

      subroutine io_dump_pl_symbap(dplfile,nbod,mass,xh,vxh,lclose,
     &                            iflgchk,rpl,rhill,j2rp2,j4rp4)

      include '../swift.inc'
      include 'io.inc'

c...    Input
      integer nbod
      real*8 mass(nbod),rpl(nbod),j2rp2,j4rp4
      real*8 xh(3,nbod),vxh(3,nbod),rhill(nbod)
      integer iflgchk
      character*(*) dplfile
      logical*2 lclose

c...   Internal
      integer j,ierr

c-----
c...  Executable code      

      call io_open(7,dplfile,'unknown','formatted',ierr)

      write(7,*) nbod

      if(btest(iflgchk,5))  then ! bit 5 is set
         write(7,123) mass(1),j2rp2,j4rp4
      else
         write(7,123) mass(1)
      endif
      write(7,123) xh(1,1),xh(2,1),xh(3,1)
      write(7,123) vxh(1,1),vxh(2,1),vxh(3,1)

      do j=2,nbod
         if(lclose) then
            write(7,123) mass(j),rhill(j),rpl(j)
         else
            write(7,123) mass(j),rhill(j)
         endif
         write(7,123) xh(1,j),xh(2,j),xh(3,j)
         write(7,123) vxh(1,j),vxh(2,j),vxh(3,j)
      enddo
 123  format(3(1p1e23.16,1x))

      close(unit = 7)
      return
      end    ! io_dump_pl_symba.f
c--------------------------------------------------------------------------

