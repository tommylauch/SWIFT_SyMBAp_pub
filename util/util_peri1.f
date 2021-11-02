c*************************************************************************
c                            UTIL_PERI1.F
c*************************************************************************
c This subroutine determines whether peri has taken place.
c Only for one particle.
c
c             Input:
c                 iflg           ==>  = 0 if first step; = 1 not (int scalar)
c                 xt,yt,zt       ==>  planocantric position of tp's
c                                       (real scalar)
c                 vxt,vyt,vzt    ==>  planocantric velcocities of tp's
c                                       (real scalar)
c                 massc          ==>  mass of the central body (real scalar)
c
c             Output:
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer scalar)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real scalar)
c                 tperi          ==> set to time to next or last perihelion, 
c                                    which ever is less, if isperi=0
c                                         (real scalar)
c
c
c Remarks: Based on util_peri
c Authors:  Hal Levison 
c Date:    2/13/01
c Last revision: 8/7/01

      subroutine util_peri1(iflg,xt,yt,zt,vxt,vyt,vzt,
     &     massc,isperi,peri,tperi)

      include '../swift.inc'

c...  Inputs Only: 
      integer ntp,iflg
      real*8 xt,yt,zt,massc
      real*8 vxt,vyt,vzt

c...  Outputs:
      real*8 peri,tperi
      integer isperi

c...  Internals
      integer i,ialpha
      real*8 vdotr,a,e
      real*8 capm              ! Not used

c----
c...  Executable code 

      if(iflg.eq.0) then    ! are we just setting thing up?

         vdotr = xt*vxt + yt*vyt + zt*vzt
         if (vdotr .gt. 0.d0) then
            isperi = 1
         else 
            isperi =-1
         endif

      else

         vdotr = xt*vxt + yt*vyt + zt*vzt
         if(isperi.eq.-1) then  ! was coming in
            
            if (vdotr .lt. 0.d0) then ! still coming in
               isperi = -1
            else                ! turned around
               isperi = 0
               call orbel_xv2aqt(xt,yt,zt,vxt,vyt,
     &              vzt,massc,ialpha,a,peri,capm,tperi)
            endif
            
         else
            
            if (vdotr .lt. 0.d0) then ! coming in
               isperi = -1
            else
               isperi = 1       ! going out
            endif
            
         endif
      
      endif

      return
      end    ! util_peri1
c------------------------------------------------------------------


