c*************************************************************************
c                        GETACCH_IR3.F
c*************************************************************************
c Calculate r^-3 for an array of particles
c             Input:
c                 nbod     ==>  number of massive bodies (int scalor)
c                istart    ==>  body to start with (int scalor)
c                 x        ==>  positions (real array)
c             Output:
c                 ir3       ==>  r^-3  (real array)
c                 ir        ==>  r^-1  (real array)
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 2/24/94

subroutine getacch_ir3(nbod,istart,x,ir3,ir)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod,istart
real(rk), intent(in)    :: x(:,:)

real(rk), intent(out)   :: ir3(:)
real(rk), intent(out)   :: ir(:)

integer(ik)             :: i
real(rk)                :: r2

c...  Executable code

do i=istart,nbod
   r2 = dot_product(x(1:3,i),x(1:3,i))
   ir(i) = 1.0_rk/sqrt(r2)
   ir3(i) = ir(i)/r2
enddo

return
end subroutine getacch_ir3
