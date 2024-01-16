module orbel_interface
implicit none

interface
   subroutine orbel_scget(angle,sx,cx)
   implicit none
   use swift_mod
   real(rk), intent(in)  :: angle
   real(rk), intent(out) :: sx,cx
   end subroutine orbel_scget
end interface

interface
   subroutine orbel_xv2aeq(x,vx,gmsum,ialpha,a,e,q)
   implicit none
   use swift_mod
   real(rk), intent(in)     :: x(:),vx(:),gmsum
   integer(ik), intent(out) :: ialpha
   real(rk), intent(out)    :: a,e,q
   end subroutine orbel_xv2aeq
end interface

interface
   subroutine orbel_xv2el(x,vx,gmsum,ialpha,a,e,inc,capom,omega,capm)
   implicit none
   use swift_mod
   real(rk), intent(in)     :: x(:),vx(:),gmsum
   integer(ik), intent(out) :: ialpha
   real(rk), intent(out)    :: a,e,inc,capom,omega,capm
   end subroutine orbel_xv2el
end interface

end module orbel_interface
