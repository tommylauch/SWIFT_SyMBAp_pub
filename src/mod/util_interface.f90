module util_interface
implicit none

interface
   subroutine util_exit(iflg)
   use swift_mod
   implicit none
   integer, intent(in) :: iflg
   end subroutine util_exit
end interface

interface
   subroutine util_hills1(msun,mpl,xh,vxh,rhill) 
   use swift_mod
   implicit none
   real(rk), intent(in)       :: msun,mpl,xh(:),vxh(:)
   real(rk), intent(out)      :: rhill 
   end subroutine util_hills1
end interface

interface
   subroutine util_hills(nbod,mass,xh,vxh,r2hill)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: r2hill(:)
   end subroutine util_hills
end interface

interface
   subroutine util_mass_peri(iflg,nbod,x,vx,mass,isperi,peri,lperi)
   use swift_mod
   implicit none
   integer(ik), intent(in)  :: nbod,iflg
   real(rk), intent(in)     :: mass(:),x(:,:),vx(:,:)
   real(rk), intent(out)    :: peri(:)
   integer(ik), intent(out) :: isperi(:)
   logical(ik), intent(out) :: lperi(:)
   end subroutine util_mass_peri
end interface

end module util_interface
