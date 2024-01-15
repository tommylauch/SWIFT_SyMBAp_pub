module swift_mod
implicit none                                                           ! you got it baby

use, intrinsic :: iso_fortran_env, only: ik => integer(ik), parameter  ::32, rk => real64

!...    Version of Swift
real(rk), parameter     :: VER_NUM = 2.0d0

!...    Maximum array size, note your stack size!
integer(ik), parameter  ::  NPLMAX = 10001                              ! max number of particles
integer(ik), parameter  ::  NTPMAX = 10001                              ! max number of test particles

!...    Size of the test particle integer status flag
integer(ik), parameter  :: NSTATP = 3                                   ! Number of status parameters
integer(ik), parameter  :: NSTAT = NSTATP + NPLMAX - 1                  ! Number of status parameters

!...    Size of the test particle integer status flag
integer(ik), parameter  :: NSTATR = NSTAT                               ! io_init_tp assumes NSTAT==NSTATR

!...    convergence criteria for danby
real(rk), parameter     :: DANBYAC = 1.0e-14_rk, DANBYB = 1.0e-13_rk

!...    loop limits in the Laguerre attempts
integer(ik), parameter  :: NLAG1 = 50, NLAG2 = 400

!...    A small number
real(rk), parameter     :: TINY = 4.0e-15_rk

!...    trig stuff
real(rk), parameter     :: PI = 3.14159265358979324_rk
real(rk), parameter     :: TWOPI = 2.0_rk*PI
real(rk), parameter     :: TWOPISQ = TWOPI**2
real(rk), parameter     :: PIBY2 = 0.5_rk*PI
real(rk), parameter     :: PI3BY2 = 1.5_rk*PI
real(rk), parameter     :: DEGRAD = 180.0_rk/PI
real(rk), parameter     :: ONETHRD = 1.0_rk/3.0_rk

end module swift_mod
