c*************************************************************************
c                        SYMBA5.INC
c*************************************************************************
C Include file for the SyMBA subroutines
c
c Remarks:  copied from symba.inc
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 

c...  Maximum number of encounters
	integer NENMAX
c	parameter(NENMAX=NTPMAX*100)   ! Must be less than 32767
	parameter(NENMAX=32767)   ! Must be less than 32767

c...	scale factor for hill's sphere to take shorter time step
       real*8 RHSCALE
       parameter (RHSCALE=6.5)

c...   Ratio of shell radii squared
	real*8 RSHELL
	parameter (RSHELL = 0.48075d0)   ! should be someting faster than RSHELL ~ NTENC^(-2/3) 

c..    ratio of the number of time steps in the adjoining shells 
       integer NTENC
       parameter (NTENC=4)


