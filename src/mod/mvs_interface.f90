module mvs_interface
implicit none

interface
   subroutine drift_dan(mu,x0,vx0,dt0,iflg)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: mu,dt0
   real(rk), intent(inout)  :: x0(:),vx0(:)
   integer(ik), intent(out) :: iflg
   end subroutine drift_dan
end interface

interface
   subroutine drift_kepmd(dm,es,ec,x,s,c)
   use swift_mod
   implicit none
   real(rk), intent(in)  :: dm,es,ec
   real(rk), intent(out) :: x,s,c
   end subroutine drift_kepmd
end interface

interface
   subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: dt,r0,mu,alpha,u
   real(rk), intent(out)    :: fp,c1,c2,c3
   integer(ik), intent(out) :: iflg
   end subroutine drift_kepu
end interface

interface
   subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: dt,r0,mu,alpha,u,s
   real(rk), intent(out)    :: f
   end subroutine drift_kepu_fchk
end interface

interface
   subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)
   use swift_mod
   implicit none
   real(rk), intent(inout) :: x
   real(rk), intent(out)   :: c0,c1,c2,c3
   end subroutine drift_kepu_stumpff
end interface

interface
   subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)
   use swift_mod
   implicit none
   real(rk), intent(in)    :: dt,r0,mu,alpha,u
   real(rk), intent(inout) :: s
   end subroutine drift_kepu_guess
end interface

interface
   subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: dt,r0,mu,alpha,u
   real(rk), intent(inout)  :: s
   real(rk), intent(out)    :: fp,c1,c2,c3
   integer(ik), intent(out) :: iflg
   end subroutine drift_kepu_lag
end interface

interface
   subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: dt,r0,mu,alpha,u
   real(rk), intent(inout)  :: s
   real(rk), intent(out)    :: fp,c1,c2,c3
   integer(ik), intent(out) :: iflgn
   end subroutine drift_kepu_new
end interface

interface
   subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: dt,r0,mu,alpha,u
   integer(ik), intent(out) :: iflg
   real(rk), intent(out)    :: s
   end subroutine drift_kepu_p3solve
end interface

interface
   subroutine drift_one(mu,x,vx,dt,iflg)
   use swift_mod
   implicit none
   real(rk), intent(in)     :: mu,dt
   real(rk), intent(inout)  :: x(:),vx(:)
   integer(ik), intent(out) :: iflg
   end subroutine drift_one
end interface

interface
   subroutine getacch_ir3(nbod,istart,x,ir3,ir)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod,istart
   real(rk), intent(in)    :: x(:,:)
   real(rk), intent(out)   :: ir3(:),ir(:)
   end subroutine getacch_ir3
end interface

interface
   subroutine kickvh(nbod,vxh,axh,dt) 
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod
   real(rk), intent(in)    :: axh(:,:),dt
   real(rk), intent(inout) :: vxh(:,:)
   end subroutine kickvh
end interface

end module mvs_interface
