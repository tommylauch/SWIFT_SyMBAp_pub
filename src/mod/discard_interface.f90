module discard_interface
implicit none

interface
   subroutine discard_massive5p(time,nbod,mass,xh,vxh,rmin,rmax,rmaxu,    &
           qmin,rpl,rhill,mergelst,mergecnt,iecnt,eoff,i1st)
   use swift_mod
   implicit none
   real(rk), intent(in)       :: time
   real(rk), intent(in)       :: rmin,rmax,rmaxu,qmin
   integer(ik), intent(in)    :: mergecnt,iecnt(:)
   integer(ik), intent(inout) :: nbod,i1st,mergelst(:,:)
   real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(inout)    :: eoff,rpl(:),rhill(:)
   end subroutine discard_massive5p

   subroutine discard_mass_merge5p_mtiny(time,ip1,ip2,mass,xh,vxh,rpl, &
                                         eoff,ielc,ielst)
   use swift_mod
   implicit none
   real(rk), intent(in)       :: time
   integer(ik), intent(inout) :: ip1,ip2
   real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:),rpl(:),eoff
   integer(ik), intent(inout) :: ielst(:,:),ielc
   end subroutine discard_mass_merge5p_mtiny

   subroutine discard_mass_peri5p(time,nbod,iecnt,mass,xh,                &
                               vxh,qmin,iwhy,isperi)
   use swift_mod
   implicit none
   integer(ik), intent(in)    :: nbod,iecnt(:)
   real(rk), intent(in)       :: mass(:),time,qmin
   real(rk), intent(in)       :: xh(:,:),vxh(:,:)
   integer(ik), intent(inout) :: iwhy(:)
   integer(ik), intent(inout) :: isperi(:)
   end subroutine discard_mass_peri5p

   subroutine discard_mass_reorder5(ip,nbod,mass,xh,vxh,rpl,rhill,isperih)
   use swift_mod
   implicit none
   integer(ik), intent(in)    :: ip
   integer(ik), intent(inout) :: nbod
   real(rk), intent(inout)    :: mass(:),xh(:,:)
   real(rk), intent(inout)    :: vxh(:,:),rpl(:)
   real(rk), intent(inout)    :: rhill(:)
   integer(ik), intent(inout) :: isperih(:)
   end subroutine discard_mass_reorder5
end interface

end module discard_interface
