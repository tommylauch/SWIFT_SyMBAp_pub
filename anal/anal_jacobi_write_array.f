c*************************************************************************
c                          ANAL_JACOBI_WRITE.F
c*************************************************************************
c Writes the mean and maximum absolute value of the change in jacobi
c to the screen as a function of time.  Writes the value of the first 10 
c test particles to a file called jacobi.out (unit=iu)
c
c      Input:
c            t              ==>  current time
c            nbod           ==>  number of massive bodies (int scalar)
c            ntp            ==>  number of tp (int scalar)
c            mass           ==>  mass of bodies (real array)
c            xh,yh,zh       ==>  current position in helio coord 
c                               (real arrays)
c            vxh,vyh,vzh    ==>  current velocity in helio coord 
c                               (real arrays)
c            xht,yht,zht    ==>  current tp position in helio coord 
c                               (real arrays)
c            vxht,vyht,vzht ==>  current tp velocity in helio coord 
c                               (real arrays)
c            istat          ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c            ipl            ==>  Planet to take jacobi with respect to 
c            iu             ==>  unit to write to
c            fopenstat      ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c
c Remarks: If the particle is not active a value of -100 is written out 
c Authors:  Hal Levison 
c Date:    3/4/93
c Last revision: 10/3/96

      subroutine anal_jacobi_write_array(t,nbod,ntp,mass,xh,vxh,
     &    xht,vxht,istat,ipl,iu,fopenstat)


      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp,ipl,iu
      integer istat(ntp)
      real*8 mass(nbod),t
      real*8 xh(3,nbod),vxh(3,nbod)
      real*8 xht(3,ntp),vxht(3,ntp)
      character*80 fopenstat

c...  Internals
      integer i1st,i,icnt,nw
      real*8 xb(3,NPLMAX),vxb(3,NPLMAX)
      real*8 xbt(3,NTPMAX),vxbt(3,NTPMAX)
      real*8 jac,jac0(NTPMAX),djmean,djmax,dj(NTPMAX)
      real*8 gmsum,energy,aplh,omega,fac
      real*8 omegax(3),msys

      data i1st/0/
      save i1st,jac0

c----
c...  Executable code 

      nw = min0(ntp,10)

c...   Compute ang. mom. vector for Sun-planet relative orbit
      gmsum = mass(1) + mass(ipl)
      energy = 0.5d0*sum(vxh(:,ipl)**2)
      energy = energy - gmsum/sqrt(sum(xh(:,ipl)**2))
      aplh = -0.5d0*gmsum/energy
      omega = sqrt(gmsum/(aplh**3))
      omegax(1) = xh(2,ipl)*vxh(3,ipl) - xh(3,ipl)*vxh(2,ipl)
      omegax(2) = xh(3,ipl)*vxh(1,ipl) - xh(1,ipl)*vxh(3,ipl)
      omegax(3) = xh(1,ipl)*vxh(2,ipl) - xh(2,ipl)*vxh(1,ipl)
      fac = omega/sqrt(sum(omegax**2))
      omegax = fac*omegax

c...  put things in bary
      call coord_h2b_array(nbod,mass,xh,vxh,xb,vxb,msys)   

      call coord_h2b_tp_array(ntp,xht,vxht,xb(:,1),vxb(:,1),xbt,vxbt)

      if(i1st.eq.0) then

         do i=1,ntp
            call anal_jacobi_array(mass(1),mass(ipl),omegax,
     &      xbt(:,i),vxbt(:,i),xb(:,1),xb(:,ipl),jac0(i))
            dj(i) = 0.0
         enddo

         call io_jacobi_write(i1st,t,jac0,dj,nw,iu,fopenstat)

         i1st = 1

      else

         icnt = 0
         djmean = 0.0
         djmax = 0.0
         do i=1,ntp
            if(istat(i).eq.0) then
               call anal_jacobi_array(mass(1),mass(ipl),omegax,
     &              xbt(:,i),vxbt(:,i),xb(:,1),xb(:,ipl),jac)
               icnt = icnt + 1
               dj(i) = jac/jac0(i) - 1.d0
               djmean = djmean + abs(dj(i))
               djmax = dmax1(djmax,abs(dj(i)))
            else
               dj(i) = -100.0
            endif
         enddo
         djmean = djmean/float(icnt)
         write(*,1) djmean,djmax
 1       format(5x,'mean |dj/j|, max |dj/j|,',2(2x,1p1e12.5))

         call io_jacobi_write(i1st,t,jac0,dj,nw,iu,fopenstat)
         
      endif

      return
      end                       ! anal_jacobi_write
c-------------------------------------------------------------------------
