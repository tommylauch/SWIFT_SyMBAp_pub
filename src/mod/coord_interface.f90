module coord_interface
implicit none

interface
   subroutine coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: msys,xb(:,:),vxb(:,:)
   end subroutine coord_h2b
end interface

interface
   subroutine coord_vb2h(nbod,mass,vxb,vxh)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:)
   real(rk), intent(inout) :: vxb(:,:)
   real(rk), intent(out)   :: vxh(:,:)
   end subroutine coord_vb2h
end interface

interface
   subroutine coord_vh2b(nbod,mass,vxh,vxb)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),vxh(:,:)
   real(rk), intent(out)   :: vxb(:,:)
   end subroutine coord_vh2b
end interface

end module coord_interface
