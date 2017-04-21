function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,flux,   &
 ferr,exptime,ntt,tobs,omc,sptmodel,nwvc)
use precision
implicit none
!import vars
integer :: nwv,nobs,nplanet,npars,nwvc
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: time,flux,ferr,exptime,tobs,omc,sptmodel
!function return type
real(double) :: loglikelihood
!local vars
integer :: i,j,nld1,nld(4)
real(double) Pi,tPi,ll1,ll2,ll3,ld(4),prior
real(double),allocatable, dimension(:,:) :: ferr2,smf,sptmodel1

interface
   subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,  &
      exptime,ntt,tobs,omc,sptmodel,nwvc)
      use precision
      implicit none
      integer :: nplanet,npars,nwv,nobs,nwvc
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: sptmodel,tobs,omc,time,  &
       exptime
   end subroutine sptransitmodel
end interface

!Physical Constants
Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi

500 format(10000(1PE17.10,1X))

!make a transit-model to compare to the data
allocate(sptmodel1(nwv,nobs))
call sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,   &
   ntt,tobs,omc,sptmodel1,nwvc)
if(nwvc.eq.0)then !entire model was calculated
   sptmodel=sptmodel1
else  !only a single bandpass was calculated
   sptmodel(nwvc,:)=sptmodel1(nwvc,:)
endif

!begin calculating log loglikelihood
!this assumes all data is valid.
ll1=0.0d0!npars*nwv*log(tPi) !number of observations

allocate(ferr2(nwv,nobs))
ferr2=ferr*ferr
ll2=0.0d0!Sum(log(ferr2)) !vectorized (seems to be faster)

allocate(smf(nwv,nobs))
smf=sptmodel-flux
ll3=Sum( smf*smf/(ferr2)) !vectorized

loglikelihood=-0.5*(ll1+ll2+ll3)

!need to add Priors
prior=1.0d0
!solerr()

!apply prior
loglikelihood=loglikelihood*prior

return
end function loglikelihood
