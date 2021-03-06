c*************************************************************************
c                            DISCARD_MASS_MERGE5P_MTINY.F
c*************************************************************************
c Merge two massive bodies
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 ip1,ip2       ==>  planets to merge (real scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 eoff          ==> Amount of energy lost due to discards
c                                          (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c             Output:
c                 mass          ==>  recalculated mass of bodies (real array)
c                 xh,yh,zh      ==>  recalculated position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  recalculated pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  recalculated physical sizes of a planet.
c                                    (real array)
c                 eoff          ==> Updated amount of energy lost from discards
c                                          (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c
c Remarks: Based on discard_mass_merge5
c
c Authors:  Hal Levison 
c Date:    12/16/09
c Last revision:  01/09/22 energy offset ignored for SyMBAp

      subroutine discard_mass_merge5p_mtiny(time,nbod,nbodm,ip1,ip2,
     &           mass,xh,yh,zh,vxh,vyh,vzh,rpl,eoff,ielc,ielst,ntpmaxsq)


      include '../swift.inc'

c...  Inputs: 
      integer ip1,ip2
      real*8 time

c...  Input and Output
      integer nbod,nbodm,ntpmaxsq
      real*8 mass(NTPMAX),xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX),rpl(NTPMAX)
      real*8 eoff
      integer ielst(2,ntpmaxsq),ielc

c...  internal
      real*8 mtot,m1,r1
      real*8 x1,y1,z1
      real*8 vx1,vy1,vz1
      real*8 m2,r2
      real*8 x2,y2,z2
      real*8 vx2,vy2,vz2
      integer itmp,j,i
      real*8 j2rp2,j4rp4

c-----
c...  Executable code 

      j2rp2=0.0d0
      j4rp4=0.0d0

      if( mass(ip2).gt.mass(ip1) ) then
         itmp = ip1
         ip1 = ip2
         ip2 = itmp
      endif

      write(*,*) ' Merging particles ',ip1, ' and ', ip2,' at t= ',time

      x1 = xh(ip1)
      y1 = yh(ip1)
      z1 = zh(ip1)
      vx1 = vxh(ip1)
      vy1 = vyh(ip1)
      vz1 = vzh(ip1)
      m1 = mass(ip1)
      r1 = rpl(ip1)

      x2 = xh(ip2)
      y2 = yh(ip2)
      z2 = zh(ip2)
      vx2 = vxh(ip2)
      vy2 = vyh(ip2)
      vz2 = vzh(ip2)
      m2 = mass(ip2)
      r2 = rpl(ip2)

c... Note:  I am just putting these guys together here, which is
c...        clearly wrong.  I should integrate back to the time
c...        of close approach.

      mtot = mass(ip1)+mass(ip2)

      rpl(ip1) = ( r1**3 + r2**3 )**(1.0d0/3.0d0)
      vxh(ip1) = (mass(ip1)*vx1 + mass(ip2)*vx2)/mtot
      vyh(ip1) = (mass(ip1)*vy1 + mass(ip2)*vy2)/mtot
      vzh(ip1) = (mass(ip1)*vz1 + mass(ip2)*vz2)/mtot
      mass(ip1) = mtot

c..   Put in zeros for the rest the second particle
      xh(ip2) = xh(ip2)*1.0d10   ! so danby does not fail
      yh(ip2) = yh(ip2)*1.0d10 
      zh(ip2) = zh(ip2)*1.0d10 
      vxh(ip2) = 0.0d0
      vyh(ip2) = 0.0d0
      vzh(ip2) = 0.0d0
      mass(ip2) = 0.0d0
      rpl(ip2) = 0.0d0

c..   Remove any encounters with ip2
      j = 1
      do while(j.le.ielc)
         if( (ielst(1,j).eq.ip2) .or. (ielst(2,j).eq.ip2) ) then
            do i=j+1,ielc
               ielst(1,i-1) = ielst(1,i)
               ielst(2,i-1) = ielst(2,i)
            enddo
            ielc = ielc - 1
         else
            j = j + 1
         endif
      enddo

      eoff = 0.d0

      call io_discard_merge(time,ip1,ip2,m1,r1,x1,y1,z1,vx1,vy1,
     &     vz1,m2,r2,x2,y2,z2,vx2,vy2,vz2,
     &     mass(ip1),rpl(ip1),xh(ip1),yh(ip1),zh(ip1),
     &     vxh(ip1),vyh(ip1),vzh(ip1))

      return
      end
