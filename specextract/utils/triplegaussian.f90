function triplegaussian(nfit,sol,x)
!returns the sum of three Gaussians and an offset
use precision
implicit none
integer :: nfit,i,np,nline
real(double) :: triplegaussian,x
real(double), dimension(nfit) :: sol

triplegaussian=sol(1)
do i=1,3
   np=3*(i-1)
   triplegaussian=triplegaussian+sol(np+2)*exp(-(x-sol(np+3))**2.0d0/   &
      (2.0d0*sol(np+4)*sol(np+4)))
enddo

return
end function triplegaussian
