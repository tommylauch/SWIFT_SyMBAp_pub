module coord_interface
implicit none

interface
   subroutine coord_h2b(nbod,mass,xh,vxh,xb,vxb,msys)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: msys,xb(:,:),vxb(:,:)
   end subroutine coord_h2b
end interface

interface
   subroutine coord_vb2h(nbod,mass,vxb,vxh)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),vxb(:,:)
   real(rk), intent(out)   :: vxh(:,:)
   end subroutine coord_vb2h
end interface

interface
   subroutine coord_vh2b(nbod,mass,vxh,vxb)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),vxh(:,:)
   real(rk), intent(out)   :: vxb(:,:)
   end subroutine coord_vh2b
end interface

end module coord_interface
