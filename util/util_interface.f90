module util_interface
implicit none

interface
   subroutine util_exit(iflg)
   implicit none
   use swift_mod
   integer, intent(in) :: iflg
   end subroutine util_exit
end interface

interface
   subroutine util_hills1(msun,mpl,xh,vxh,rhill) 
   implicit none
   use swift_mod
   real(rk), intent(in)       :: msun,mpl,xh(:),vxh(:)
   real(rk), intent(out)      :: rhill 
   end subroutine util_hills1
end interface

interface
   subroutine util_hills(nbod,mass,xh,vxh,r2hill)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: r2hill(:)
   end subroutine util_hills
end interface

interface
   subroutine util_mass_peri(iflg,nbod,x,vx,mass,isperi,peri,lperi)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: nbod,iflg
   real(rk), intent(in)     :: mass(:),x(:,:),vx(:,:),gm
   real(rk), intent(out)    :: peri(:)
   integer(ik), intent(out) :: isperi(:)
   logical(ik), intent(out) :: lperi(:)
   end subroutine util_mass_peri
end interface

end module util_interface
