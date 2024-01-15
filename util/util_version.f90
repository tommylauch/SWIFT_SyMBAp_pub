subroutine util_version

use swift_mod

c-----
c...  Executable code 

write(*,*)            '-----------------------------------------------'
write(*,'(a,f3.1,a)') '------------- SyMBAp: version ',VER_NUM,' -------------'
write(*,*)            '-----------------------------------------------'

return
end subroutine util_version
