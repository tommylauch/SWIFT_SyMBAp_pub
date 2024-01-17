module helio_interface
implicit none

interface
   subroutine helio_drift(nbod,mass,xh,vxb,dt)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),dt
   real(rk), intent(inout) :: xh(:,:),vxb(:,:)
   end subroutine helio_drift

   subroutine helio_lindrift(nbod,mass,vxb,dt,xh,ptx)
   use swift_mod
   implicit none
   integer(ik), intent(in)  :: nbod
   real(rk), intent(in)     :: mass(:),dt,vxb(:,:)
   real(rk), intent(inout)  :: xh(:,:)
   real(rk), intent(out)    :: ptx(:)
   end subroutine helio_lindrift
end interface

end module helio_interface
