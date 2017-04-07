subroutine EstZpt(npars,sol,solerr,nwv,nobs,time,flux,exptime,ntt,tobs, &
 omc)
use precision
implicit none
integer :: npars,nwv,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol,time,exptime
real(double), dimension(:,:) :: solerr,flux,tobs,omc

return
end program estzpt
