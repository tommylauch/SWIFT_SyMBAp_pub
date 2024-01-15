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
c                 xh            ==>   position in helio coord 
c                                    (real arrays)
c                 vxh           ==>   pl vel in helio coord 
c                                    (real arrays)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 eoff          ==> Amount of energy lost due to discards
c                                          (real scalar)
c                ielc           ==>  number of encounters (integer scalar)
c                ielst          ==>  list of ecnounters (2D integer array)
c             Output:
c                 mass          ==>  recalculated mass of bodies (real array)
c                 xh            ==>  recalculated position in helio coord 
c                                    (real arrays)
c                 vxh           ==>  recalculated pl vel in helio coord 
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
     &           mass,xh,vxh,rpl,eoff,ielc,ielst,ntpmaxsq)


      include '../swift.inc'

c...  Inputs: 
      integer ip1,ip2
      real*8 time

c...  Input and Output
      integer nbod,nbodm,ntpmaxsq
      real*8 mass(NTPMAX),xh(3,NTPMAX),vxh(3,NTPMAX),rpl(NTPMAX),eoff
      integer ielst(2,ntpmaxsq),ielc

c...  internal
      real*8 mtot,m(2),r(2),x(3,2),vx(3,2)
      integer itmp,j,i
      real*8 j2rp2,j4rp4

c-----
c...  Executable code 

      j2rp2 = 0.0d0
      j4rp4 = 0.0d0

      if( mass(ip2).gt.mass(ip1) ) then
         itmp = ip1
         ip1 = ip2
         ip2 = itmp
      endif

      write(*,*) ' Merging particles ',ip1, ' and ', ip2,' at t= ',time

      x(:,1) = xh(:,ip1)
      vx(:,1) = vxh(:,ip1)
      m(1) = mass(ip1)
      r(1) = rpl(ip1)

      x(:,2) = xh(:,ip2)
      vx(:,2) = vxh(:,ip2)
      m(2) = mass(ip2)
      r(2) = rpl(ip2)

c... Note:  I am just putting these guys together here, which is
c...        clearly wrong.  I should integrate back to the time
c...        of close approach.

      mtot = mass(ip1)+mass(ip2)

      rpl(ip1) = ( r(1)**3 + r(2)**3 )**(1.0d0/3.0d0)
      vxh(:,ip1) = (mass(ip1)*vx(:,1) + mass(ip2)*vx(:,2))/mtot
      mass(ip1) = mtot

c..   Put in zeros for the rest the second particle
      xh(:,ip2) = xh(:,ip2)*1.0d10   ! so danby does not fail
      vxh(:,ip2) = 0.0d0
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

      call io_discard_merge_symbap(time,ip1,ip2,m,r,x,vx,
     &                mass(ip1),rpl(ip1),xh(:,ip1),vxh(:,ip1))

      return
      end
