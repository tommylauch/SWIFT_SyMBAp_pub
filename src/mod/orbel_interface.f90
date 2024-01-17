module orbel_interface
implicit none

interface
   subroutine orbel_scget(angle,sx,cx)
   use swift_mod
   implicit none
   real(rk), intent(in)  :: angle
   real(rk), intent(out) :: sx,cx
   end subroutine orbel_scget

   subroutine orbel_xv2aeq(x,vx,gmsum,ialpha,a,e,q)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: x(:),vx(:),gmsum
   integer(ik), intent(out) :: ialpha
   real(rk), intent(out)    :: a,e,q
   end subroutine orbel_xv2aeq

   subroutine orbel_xv2el(x,vx,gmsum,ialpha,a,e,inc,capom,omega,capm)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: x(:),vx(:),gmsum
   integer(ik), intent(out) :: ialpha
   real(rk), intent(out)    :: a,e,inc,capom,omega,capm
   end subroutine orbel_xv2el
end interface

end module orbel_interface
