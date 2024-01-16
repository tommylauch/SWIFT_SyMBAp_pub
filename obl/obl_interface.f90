module obl_interface
implicit none

interface
   subroutine obl_acc(nbod,mass,j2rp2,j4rp4,xh,irh,aoblx)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: j2rp2,j4rp4
   real(rk), intent(in)    :: mass(:),xh(:,:),irh(:)
   real(rk), intent(out)   :: aoblx(:,:)
   end subroutine obl_acc
end interface

interface
   subroutine obl_pot(nbod,mass,j2rp2,j4rp4,xh,irh,oblpot)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:)
   real(rk), intent(in)    :: j2rp2,j4rp4
   real(rk), intent(in)    :: xh(:,:),irh(:)
   real(rk), intent(out)   :: oblpot
   end subroutine obl_pot
end interface

end module obl_interface
