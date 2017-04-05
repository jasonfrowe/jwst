subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,     &
 sptmodel)
use precision
implicit none
integer :: nplanet,npars,nwv,nobs
real(double), dimension(:) :: sol,time
real(double), dimension(:,:) :: solrange,sptmodel

return
end subroutine sptransitmodel
