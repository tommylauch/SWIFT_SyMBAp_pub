module discard_interface
implicit none

interface
   subroutine discard_massive5p(time,dt,nbod,mass,xh,vxh,rmin,rmax,rmaxu, &
             qmin,lclose,rpl,rhill,isenc,mergelst,mergecnt,iecnt,eoff,i1st)
   implicit none
   use swift_mod
   real(rk), intent(in)       :: time,dt
   integer(ik), intent(in)    :: isenc
   real(rk), intent(in)       :: rmin,rmax,rmaxu,qmin
   logical(ik), intent(in)    :: lclose
   integer(ik), intent(in)    :: mergelst(:,:),mergecnt,iecnt(:)
   integer(ik), intent(inout) :: nbod,i1st
   real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:)
   real(rk), intent(inout)    :: eoff,rpl(:),rhill(:)
   end subroutine discard_massive5p
end interface

interface
   subroutine discard_mass_merge5p_mtiny(time,nbod,nbodm,ip1,ip2,         &
           mass,xh,vxh,rpl,eoff,ielc,ielst)
   implicit none
   use swift_mod
   integer(ik), intent(in)    :: ip1,ip2
   real(rk), intent(in)       :: time
   integer(ik), intent(inout) :: nbod,nbodm
   real(rk), intent(inout)    :: mass(:),xh(:,:),vxh(:,:),rpl(:),eoff
   integer(ik), intent(inout) :: ielst(:,:),ielc
   end subroutine discard_mass_merge5p_mtiny
end interface

interface
   subroutine discard_mass_peri5p(time,nbod,iecnt,mass,xh,                &
                               vxh,qmin,iwhy,isperi)
   implicit none
   use swift_mod
   integer(ik), intent(in)    :: nbod,iecnt(:)
   real(rk), intent(in)       :: mass(:),time,qmin
   real(rk), intent(in)       :: xh(:,:),vxh(:,:)
   integer(ik), intent(inout) :: iwhy(:)
   integer(ik), intent(inout) :: isperi(:)
   end subroutine discard_mass_peri5p
end interface

interface
   subroutine discard_mass_reorder5(ip,nbod,mass,xh,vxh,rpl,rhill,isperih)
   implicit none
   use swift_mod
   integer(ik), intent(in)    :: ip
   integer(ik), intent(inout) :: nbod
   real(rk), intent(inout)    :: mass(:),xh(:,:)
   real(rk), intent(inout)    :: vxh(:,:),rpl(:)
   real(rk), intent(inout)    :: rhill(:)
   integer(ik), intent(inout) :: isperih(:)
   end subroutine discard_mass_reorder5
end interface

end module discard_interface
