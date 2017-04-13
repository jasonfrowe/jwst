function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,flux,   &
 ferr,exptime,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: nwv,nobs,nplanet,npars
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: time,flux,ferr,exptime,tobs,omc
!function return type
real(double) :: loglikelihood
!local vars
integer :: i,j
real(double) Pi,tPi,ll1,ll2,ll3
real(double),allocatable, dimension(:,:) :: sptmodel

interface
   subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,  &
      exptime,ntt,tobs,omc,sptmodel)
      use precision
      implicit none
      integer :: nplanet,npars,nwv,nobs
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

!make a transit-model to compare to the data
allocate(sptmodel(nwv,nobs)) !array to hold the spectral transit model
call sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,   &
   ntt,tobs,omc,sptmodel)

!begin calculating log loglikelihood
!this assumes all data is valid.
ll1=npars*nwv*log(tPi) !number of observations

ll2=0.0d0
do i=1,nwv    !sum over uncertainties
   do j=1,nobs
      ll2=ll2+log(ferr(i,j)*ferr(i,j))
   enddo
enddo

ll3=0.0d0
do i=1,nwv
   do j=1,nobs
      ll3=ll3+(sptmodel(i,j)-flux(i,j))*(sptmodel(i,j)-flux(i,j))/      &
         (ferr(i,j)*ferr(i,j))
   enddo
enddo

loglikelihood=-0.5*(ll1+ll2+ll3)

!need to add Priors and constraints.
! * valid limb-darkening
! * 0 < e < 1

return
end function loglikelihood
