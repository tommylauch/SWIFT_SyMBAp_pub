c*************************************************************************
c                            UTIL_MASS_PERI.F
c*************************************************************************
c This subroutine determines whether peri of a planet has taken place
c
c             Input:
c                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of bodies (int scalar)
c                 x,y,z          ==>  heliocentric position of planets
c                                       (real arrays)
c                 vx,vy,vz       ==>  heliocentric velcocities of planets
c                                       (real arrays)
c                 mass           ==>  mass of the bodies (real array)
c
c             Output:
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c                lperi           ==> set to .true. if isperi=0
c                                         (logical*2 array)
c
c
c Remarks: Based on util_peri.f
c Authors:  Hal Levison 
c Date:    12/30/96
c Last revision: 

      subroutine util_mass_peri_symbap(iflg,nbod,x,vx,mass,isperi,peri,
     &                                 lperi)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,iflg
      real*8 mass(*),x(3,*),vx(3,*),gm

c...  Outputs:
      real*8 peri(*)
      integer isperi(*)
      logical*2 lperi(*)

c...  Internals
      integer i,ialpha
      real*8 vdotr,a,e

c----
c...  Executable code 

      if(iflg.eq.0) then    ! are we just setting thing up?
         do i=2,nbod
            vdotr = x(1,i)*vx(1,i)+x(2,i)*vx(2,i)+x(3,i)*vx(3,i)
            if (vdotr .gt. 0.d0) then
               isperi(i) = 1
            else 
               isperi(i) = -1
            endif
         enddo
      else
         do i=2,nbod
            vdotr = x(1,i)*vx(1,i) + x(2,i)*vx(2,i) + x(3,i)*vx(3,i)
            if(isperi(i).eq.-1) then         ! was coming in
               if (vdotr .lt. 0.d0) then     ! still coming in
                  isperi(i) = -1
               else                          ! turned around
                  isperi(i) = 0
                  lperi(i) = .true.
                  gm = mass(1) + mass(i)
                  call orbel_xv2aeq_symbap(x(1:3,i),vx(1:3,i),gm,ialpha,
     &                                     a,e,peri(i))
               endif
            else

               if (vdotr .lt. 0.d0) then     ! coming in
                  isperi(i) = -1
               else
                  isperi(i) = 1              ! going out
               endif

            endif
         enddo
      endif

      return
      end    ! util_mass_peri_symbap
c------------------------------------------------------------------


