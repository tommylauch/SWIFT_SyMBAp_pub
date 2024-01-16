module anal_interface
implicit none

interface
   subroutine anal_energy_discard5(iflg,nbod,nbodm,mass,j2rp2,j4rp4,      &
                                   xh,vxh,ke,pot,energy,eltot)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: iflg,nbod,nbodm
   real(rk), intent(in)    :: mass(:),j2rp2,j4rp4
   real(rk), intent(in)    :: xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: energy,eltot(:),ke,pot
   end subroutine anal_energy_discard5
end interface

interface
   subroutine anal_energy(nbod,mass,j2rp2,j4rp4,xh,vxh,ke,pot,energy,eltot)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: mass(:),j2rp2,j4rp4
   real(rk), intent(in)    :: xh(:,:),vxh(:,:)
   real(rk), intent(out)   :: energy,eltot(:),ke,pot
   end subroutine anal_energy
end interface

interface
   subroutine anal_energy_mtiny(nbod,nbodm,mass,j2rp2,j4rp4,xh,vxh,       &
                                ke,pot,energy,eltot)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,nbodm
   real(rk), intent(in)    :: mass(:),xh(:,:),vxh(:,:),j2rp2,j4rp4
   real(rk), intent(out)   :: energy,eltot(:),ke,pot
   end subroutine anal_energy_mtiny
end interface

interface
   subroutine anal_energy_write(t,nbod,mass,j2rp2,j4rp4,xh,vxh,           &
                                iu,fopenstat,eoff)
   implicit none
   use swift_mod
   integer(ik), intent(in)        :: nbod,iu
   real(rk), intent(in)           :: mass(:),t,j2rp2,j4rp4,eoff
   real(rk), intent(in)           :: xh(:,:),vxh(:,:)
   character(len = :), intent(in) :: fopenstat
   end subroutine anal_energy_write
end interface

end module anal_interface
