subroutine EstZpt(npars,sol,solerr,solrange,nwv,nobs,time,flux,exptime, &
 ntt,tobs,omc)
use precision
implicit none
integer :: npars,nwv,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc,      &
   solrange

!check if we have one zpt per bandpass, or just one.
if(solrange(8,2)-solrange(8,1).gt.0)then !8 corresponds to zero-point

   !calculate if data is in transit or out of transit.
   !this means calculating b and using rprs.

else !else calculate for each bandpass

endif

return
end subroutine estzpt
