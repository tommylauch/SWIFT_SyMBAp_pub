module helio_interface
implicit none

interface
   subroutine helio_drift(nbod,mass,xh,vxb,dt)	
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),dt
   real(rk), intent(inout) :: xh(:,:),vxb(:,:)
   end subroutine helio_drift
end interface

interface
   subroutine helio_lindrift(nbod,mass,vxb,dt,xh,ptx)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: nbod
   real(rk), intent(in)     :: mass(:),dt,vxb(:,:)
   real(rk), intent(inout)  :: xh(:,:)
   real(rk), intent(out)    :: ptx(:)
   end subroutine helio_lindrift
end interface

end module helio_interface
