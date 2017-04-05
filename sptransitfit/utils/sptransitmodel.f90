subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,     &
 exptime,sptmodel)
use precision
implicit none
!import vars
integer :: nplanet,npars,nwv,nobs
real(double), dimension(:) :: sol,time,exptime
real(double), dimension(:,:) :: solrange,sptmodel
!local vars
integer :: nintg
real(double) :: Pi,tPi,pid2,G,Cs,fDB

!Model parameters
nintg=41 !number of samples to convolve integration time

!Physical Constants
Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi
pid2=Pi/2.0d0  !and Pi/2
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Cs=2.99792458e8 !Speed of light
fDB=1.0 !Doppler Boosting factor


return
end subroutine sptransitmodel
