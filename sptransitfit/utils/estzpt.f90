subroutine EstZpt(npars,sol,solerr,solrange,nwv,nobs,time,flux,exptime, &
 ntt,tobs,omc)
use precision
implicit none
integer :: npars,nwv,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc,solrange

!if only one zpt then...


!calculate if data is in transit or out of transit.
!this means calculating b and using rprs.



!else calculate for each bandpass

!end

return
end subroutine estzpt
