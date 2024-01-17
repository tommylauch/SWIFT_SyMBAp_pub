module io_interface
implicit none

interface
   subroutine io_discard_cleanup(t0)
   use swift_mod
   implicit none
   real(rk), intent(in) :: t0
   end subroutine io_discard_cleanup
end interface

interface
   subroutine io_discard_mass(init,time,id,m1,r1,x1,vx1,iu,iwhy,fopenstat)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: iwhy,iu,init,id
   real(rk), intent(in)           :: time,m1,r1
   real(rk), intent(in)           :: x1(:),vx1(:)
   character(len=50), intent(in)  :: fopenstat
   end subroutine io_discard_mass
end interface

interface
   subroutine io_discard_merge(time,ip1,ip2,m,r,x,vx,mn,rn,xn,vxn)
   use swift_mod
   implicit none
   integer(ik), intent(in) :: ip1,ip2
   real(rk), intent(in)    :: time
   real(rk), intent(in)    :: m(:),r(:),x(:,:),vx(:,:)
   real(rk), intent(in)    :: mn,rn,xn(:),vxn(:)
   end subroutine io_discard_merge
end interface

interface
   subroutine io_dump_param(dparfile,t,tstop,dt,dtout,dtdump,             &
                         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
   use swift_mod
   implicit none
   real(rk), intent(in)           :: t,tstop,dt
   integer(ik), intent(in)        :: iflgchk
   real(rk), intent(in)           :: dtout,dtdump
   real(rk), intent(in)           :: rmin,rmax,rmaxu,qmin
   logical(ik), intent(in)        :: lclose
   character(len=50), intent(in)  :: outfile,dparfile
   end subroutine io_dump_param
end interface

interface
   subroutine io_dump_pl(dplfile,nbod,mass,xh,vxh,lclose,                 &
                      iflgchk,rpl,rhill,j2rp2,j4rp4)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: nbod,iflgchk
   real(rk), intent(in)           :: mass(:),rpl(:),j2rp2,j4rp4
   real(rk), intent(in)           :: xh(:,:),vxh(:,:),rhill(:)
   character(len=50), intent(in)  :: dplfile
   logical(ik), intent(in)        :: lclose
   end subroutine io_dump_pl
end interface

interface
   subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: iu,i1st
   real(rk), intent(in)           :: t,energy,eltot(:)
   character(len=50), intent(in)  :: fopenstat
   end subroutine io_energy_write
end interface

interface
   subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,iflgchk,      &
                         rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)
   use swift_mod
   implicit none
   character(len=50), intent(in)   :: infile
   integer(ik), intent(out)        :: iflgchk
   real(rk), intent(out)           :: t0,tstop,dt
   real(rk), intent(out)           :: dtout,dtdump
   real(rk), intent(out)           :: rmin,rmax,rmaxu,qmin
   logical(ik), intent(out)        :: lclose
   character(len=50), intent(out)  :: outfile,fopenstat
   end subroutine io_init_param
end interface

interface
   subroutine io_init_pl(infile,lclose,iflgchk,nbod,mass,xh,vxh,rpl,rhill,&
                      j2rp2,j4rp4)
   use swift_mod
   implicit none
   character(len=50), intent(in)  :: infile
   integer(ik), intent(in)        :: iflgchk
   logical(ik), intent(in)        :: lclose
   real(rk), intent(out)          :: mass(:),rpl(:),j2rp2,j4rp4
   real(rk), intent(out)          :: xh(:,:),vxh(:,:),rhill(:)
   integer(ik), intent(out)       :: nbod
   end subroutine io_init_pl
end interface

interface
   subroutine io_open(iu,fname,fopenstat,format,ierr)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: iu
   character(len=50), intent(in)  :: fname,fopenstat,format
   integer(ik), intent(out)       :: ierr
   end subroutine io_open
end interface

interface
   subroutine io_splitname(oname,dirname,ldir,filename,lfile)
   use swift_mod
   implicit none
   character(len=50), intent(in)  :: oname
   integer(ik), intent(out)       :: ldir,lfile
   character(len=50), intent(out) :: dirname,filename
   end subroutine io_splitname
end interface

interface
   subroutine io_write_frame(time,nbod,ntp,mass,xh,vxh,xht,vxht,istat,    &
                          oname,iu,fopenstat)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: nbod,ntp,iu
   real(rk), intent(in)           :: mass(:),time
   integer(ik), intent(in)        :: istat(:,:)
   real(rk), intent(in)           :: xh(:,:),vxh(:,:)
   real(rk), intent(in)           :: xht(:,:),vxht(:,:)
   character(len=50), intent(in)  :: oname,fopenstat
   end subroutine io_write_frame
end interface

interface
   subroutine io_write_hdr(iu,time,nbod,ntp,istat) 
   use swift_mod
   implicit none
   integer(ik), intent(in) :: nbod,ntp,istat(:,:),iu
   real(rk), intent(in)    :: time
   end subroutine io_write_hdr
end interface

interface
   subroutine io_write_line(iu,id,a,e,inc,capom,omega,capm) 
   use swift_mod
   implicit none
   integer(ik), intent(in) :: iu,id
   real(rk), intent(in)    :: a,e,inc,capom,omega,capm
   end subroutine io_write_line
end interface

interface
   subroutine io_write_mass(time,nbod,mass,oname,iu,fopenstat)
   use swift_mod
   implicit none
   integer(ik), intent(in)        :: nbod,iu
   real(rk), intent(in)           :: mass(:),time
   character(len=50), intent(in)  :: oname,fopenstat 
   end subroutine io_write_mass
end interface

end module io_interface
