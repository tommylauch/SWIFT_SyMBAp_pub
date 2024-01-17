module symba5p_mod
use swift_mod
implicit none
!...  Maximum number of encounters
integer(ik), parameter :: NENMAX = 5000_rk
!...	scale factor for hill's sphere to take shorter time step
real(rk), parameter    :: RHSCALE = 6.5_rk

!...   Ratio of shell radii squared
real(rk), parameter    :: RSHELL = 0.480749856769136133_rk   ! should be someting faster than RSHELL ~ NTENC^(-2/3) 

!..    ratio of the number of time steps in the adjoining shells 
integer(ik), parameter :: NTENC = 3_ik

!..    Maximum no. of pairs in a group
integer(ik), parameter :: GRPMAX = 2000_ik
!..    Maximum no. of group
integer(ik), parameter :: GRPNMAX = NENMAX

end module symba5p_mod
