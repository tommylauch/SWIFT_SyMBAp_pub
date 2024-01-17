module rmvs_interface
implicit none

interface
   subroutine rmvs_chk_ind(xr,vxr,dt,r2crit,r2critp,iflag)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: xr(:),vxr(:),dt,r2crit,r2critp
   integer(ik), intent(out) :: iflag
   end subroutine rmvs_chk_ind
end interface

end module rmvs_interface
