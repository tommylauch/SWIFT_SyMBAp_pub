module sybam5_mod
implicit none
use swift_mod
c...  Maximum number of encounters
integer(ik), parameter :: NENMAX = 262144
c...	scale factor for hill's sphere to take shorter time step
real(rk), parameter    :: RHSCALE = 6.5_rk

c...   Ratio of shell radii squared
real(rk), parameter    :: RSHELL = 0.480749856769136133_rk   ! should be someting faster than RSHELL ~ NTENC^(-2/3) 

c..    ratio of the number of time steps in the adjoining shells 
integer(ik), parameter :: NTENC = 3_ik

end module sybam5_mod
