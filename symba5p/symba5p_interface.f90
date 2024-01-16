module symba5p_interface
implicit none

interface
   subroutine symba5p_chk(rhill,nbod,ip1,ip2,mass,xh,vxh,dt,              &
                          irec,icflg,svdotr)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: nbod,irec,ip1,ip2
   real(rk), intent(in)     :: mass(:),xh(:,:),dt
   real(rk), intent(in)     :: vxh(:,:),rhill(:)
   integer(ik), intent(out) :: icflg
   logical(ik), intent(out) :: svdotr
   end subroutine 
end interface symba5p_chk

interface
   subroutine symba5p_getacch(nbod,nbodm,mass,j2rp2,j4rp4,xh,axh,         &
                              mtiny,ielc,ielst)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,nbodm,ielst(:,:),ielc
   real(rk), intent(in)    :: mass(:),xh(:,:),j2rp2,j4rp4,mtiny
   real(rk), intent(out)   :: axh(:,:)
   end subroutine symba5p_getacch
end interface

interface
   subroutine symba5p_group(ielst,ielc,i_ie,j_ie,grpie,grppc,grpc)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: ielst(:,:),ielc
   integer, intent(inout)   :: grpie(:,:),grppc(:),grpc
   end subroutine symba5p_group
end interface

interface
   subroutine symba5p_helio_drift(nbod,ielev,irec,mass,xh,vxb,dt)	
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,irec,ielev(:)
   real(rk), intent(in)    :: mass(:),dt
   real(rk), intent(inout) :: xh(:,:),vxb(:,:)
   end subroutine symba5p_helio_drift
end interface

interface
   subroutine symba5p_helio_drift_g(nbod,ielev,irec,mass,xh,vxb,          &
                                    dt,ielc,ielst)	
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,irec
   real(rk), intent(in)    :: mass(:),dt
   integer(ik), intent(in) :: ielev(:),ielst(:,:),ielc
   real(rk), intent(inout) :: xh(:,:),vxb(:,:)
   end subroutine symba5p_helio_drift_g
end interface

interface
   subroutine symba5p_helio_getacch(iflg,nbod,nbodm,mass,j2rp2,j4rp4,     &
                                    xh,axh)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,nbodm,iflg
   real(rk), intent(in)    :: mass(:),xh(:,:),j2rp2,j4rp4
   real(rk), intent(out)   :: axh(:,:)
   end subroutine symba5p_helio_getacch
end interface

interface
   subroutine symba5p_kick(nbod,mass,irec,iecnt,ielev,                    &
                           rhill,xh,vxb,dt,sgn,ielc,ielst)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,irec
   real(rk), intent(in)    :: mass(:),dt,rhill(:),sgn
   integer(ik), intent(in) :: iecnt(:),ielev(:),ielst(:,:),ielc
   real(rk), intent(in)    :: xh(:,:)
   real(rk), intent(inout) :: vxb(:,:)
   end subroutine symba5p_kick
end interface

interface
   subroutine symba5p_merge(t,dt,nbod,nbodm,ip1,ip2,mass,xh,vxb,          &
         ireci,irecl,svdotrold,iecnt,rpl,mergelst,mergecnt,rhill,         &
         eoff,ielc,ielst)
   implicit none
   use swift_mod
   integer(ik), intent(in)    :: nbod,nbodm,ireci,irecl,ip1,ip2
   real(rk), intent(in)       :: t,dt
   logical(ik), intent(in)    :: svdotrold
   real(rk), intent(inout)    :: mass(:),xh(:,:),vxb(:,:),eoff
   real(rk), intent(inout)    :: rpl(:),rhill(:)
   integer(ik), intent(inout) :: iecnt(:),ielst(:,:),ielc
   integer(ik), intent(inout) :: mergelst(:,:),mergecnt
   integer(ik), intent(inout) :: ip1l,ip2l
   end subroutine symba5p_merge
end interface

interface
   subroutine symba5p_nbodm(nbod,mass,mtiny,nbodm)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: nbod
   real(rk), intent(in)     :: mass(nbod),mtiny
   integer(ik), intent(out) :: nbodm
   end subroutine symba5p_nbodm
end interface

interface
   subroutine symba5p_step_helio(i1st,nbod,nbodm,mass,j2rp2,j4rp4,        &
                                 xh,vxh,dt)
   implicit none
   use swift_mod
   integer(ik), intent(in) :: nbod,i1st,nbodm
   real(rk), intent(in)    :: mass(:),dt,j2rp2,j4rp4
   real(rk), intent(inout) :: xh(:,:),vxh(:,:)
   end subroutine symba5p_step_helio
end interface

interface
   subroutine symba5p_step_interp(time,iecnt,ielev,nbod,nbodm,mass,       &
         rhill,j2rp2,j4rp4,lclose,rpl,xh,vxh,dt,mergelst,mergecnt,        &
         eoff,ielc,ielst,mtiny,grpie,grppc,grpc)
   implicit none
   use swift_mod
   real(rk), intent(in)       :: mass(:),dt,j2rp2,j4rp4,time,mtiny
   integer(ik), intent(in)    :: iecnt(:),ielev(:),ielst(:,:),ielc
   integer(ik), intent(in)    :: grpie(:,:),grppc(:),grpc
   logical(ik), intent(in)    :: lclose
   integer(ik), intent(inout) :: nbod,nbodm
   real(rk), intent(inout)    :: xh(:,:),vxh(:,:),rpl(:),eoff,rhill(:)
   integer(ik), intent(out)   :: mergelst(:,:),mergecnt
   end subroutine symba5p_step_interp
end interface

interface
   subroutine symba5p_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,j4rp4,      &
                              xh,vxh,dt,lclose,rpl,isenc,                 &
                              mergelst,mergecnt,iecnt,eoff,rhill,mtiny)
   implicit none
   use swift_mod
   integer(ik), intent(in)  :: nbod,i1st,nbodm
   real(rk), intent(in)     :: mass(:),dt,time,j2rp2,j4rp4,mtiny
   logical(ik), intent(in)  :: lclose
   real(rk), intent(inout)  :: xh(:,:),vxh(:,:)
   real(rk), intent(inout)  :: rpl(:),eoff,rhill(:)
   integer(ik), intent(out) :: isenc
   integer(ik), intent(out) :: iecnt(:),ielev(:)
   integer(ik), intent(out) :: mergelst(2,:),mergecnt
   end subroutine symba5p_step_pl
end interface

interface
   recursive subroutine symba5p_step_recur(t,nbod,nbodm,mass,ireci,ilevl, &
                        iecnt,ielev,rhill,xh,vxb,lclose,rpl,              &
                        mergelst,mergecnt,dt0,eoff,svdotr,ielc,ielst)
   implicit none
   use swift_mod
   integer(ik), intent(in)    :: nbod,ireci,nbodm
   real(rk), intent(in)       :: mass(:),dt0,rhill(:),t
   integer(ik), intent(in)    :: iecnt(:),ielev(:),ielst(:,:),ielc
   logical(ik), intent(in)    :: lclose
   integer(ik), intent(inout) :: ilevl(:)
   real(rk), intent(inout)    :: xh(:,:),vxb(:,:),eoff,rpl(:)
   integer(ik), intent(inout) :: mergelst(:,:),mergecnt
   logical(ik), intent(inout) :: svdotr(:)
   end subroutine symba5p_step_recur
end interface

end module symba5p_interface
