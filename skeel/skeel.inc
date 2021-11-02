c*************************************************************************
c                        SKEEL.INC
c*************************************************************************
C Include file for the Skeel subroutines
c
c Author:  Hal Levison  
c Date:    9/23/96
c Last revision: 

c...	scale factor for hill's sphere to take shorter time step
       real*8 RHSCALE
c       parameter (RHSCALE=3.5)
       parameter (RHSCALE=6.5)

c..    ratio of the number of time steps in the adjoining shells 
       integer NTENC
       parameter (NTENC=3)

c...   Ratio of shell radii squared
	real*8 RSHELL
	parameter (RSHELL = 0.48075d0)   ! RSHELL ~ NTENC^(-2/3) 
c	parameter (RSHELL = 0.57735d0)   ! RSHELL ~ NTENC^(-1/2) 

c...   This is for the peri stuff.  We consider a 2-body orbit good if the
c...   closest peri passage is within FACQDT*dt
	real*8 FACQDT
	parameter (FACQDT = 2.0d0)  

c...   Flag that communicates that the peri fit was not good
c...   We check that this is .gt. 0.  So this must be less than 0.
	real*8 QFLG
	parameter (QFLG = -1.0d0)  
