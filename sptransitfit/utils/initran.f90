subroutine initran(seed)
use precision
implicit none
!import vars
integer :: seed
!local vars
integer :: now(3)
real(double) ran2,dumr

!use current time to generate seed
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed) !initialize random number

return
end