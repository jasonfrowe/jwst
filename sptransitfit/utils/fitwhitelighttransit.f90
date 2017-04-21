subroutine fitwhitelighttransit(npars,nplanet,sol,solerr,solrange,nwv,  &
 nobs,time,flux,ferr,exptime,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: npars,nwv,nobs,nplanet
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc
!local vars
integer :: i,ii,j,npars_wl,nwv_wl,np,nwvc,ntest
integer, allocatable, dimension(:,:) :: solrange_wl
real(double) :: dnobs
real(double), allocatable, dimension(:) :: tn,sol_wl
real(double), allocatable, dimension(:,:) :: time_wl,flux_wl,ferr_wl,   &
 exptime_wl,solerr_wl,sptmodel

!interfaces to allow assumed arrays
interface
   subroutine fittransitmodel8(npars,nplanet,sol,solerr,solrange,nwv, &
    nobs,time,flux,ferr,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: npars,nwv,nobs,nplanet
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,    &
       tobs,omc
   end subroutine fittransitmodel8
   subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,  &
      exptime,ntt,tobs,omc,sptmodel,nwvc)
      use precision
      implicit none
      integer :: nplanet,npars,nwv,nobs,nwvc
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: sptmodel,tobs,omc,time,exptime
   end subroutine sptransitmodel
end interface

!Create whitelight photometry

!allocate array space for white light photometry
allocate(time_wl(1,nobs),flux_wl(1,nobs),ferr_wl(1,nobs),               &
 exptime_wl(1,nobs),tn(nobs))

!create weighted averages across all bandpasses
dnobs=dble(nobs) !converting int to dble
do i=1,nobs !loop over time steps
   tn(i)=sum(1/ferr(:,i)) !sum of weights for weighted average
   time_wl(1,i)=Sum(time(:,i)/ferr(:,i)) !sum of weighted times
   flux_wl(1,i)=Sum(flux(:,i)/ferr(:,i)) !sum of weighted fluxes
   ferr_wl(1,i)=dnobs  !number of observations
   exptime_wl(1,i)=Sum(exptime(:,i)/ferr(:,i))
enddo
time_wl(1,:)=time_wl(1,:)/tn !weighted time-stamp average
flux_wl(1,:)=flux_wl(1,:)/tn !weighted flux average
ferr_wl(1,:)=(ferr_wl(1,:)**0.5)/tn !propagated uncertainty on mean
exptime_wl(1,:)=exptime_wl(1,:)/tn !weighted exptime average

!create input model solution guess for whitelight model
!make sol,solerr,solrange
npars_wl=8+10*nplanet !number of parameters for whitelight model
allocate(sol_wl(npars_wl),solerr_wl(npars_wl,4),solrange_wl(npars_wl,2))

do i=1,8 !loop over global parameters
   sol_wl(i)=sol(solrange(i,1))
   solerr_wl(i,:)=solerr(solrange(i,1),:)
   solrange_wl(i,:)=i
enddo

!planet parameters
do np=1,nplanet !loop over each planet in transit model
   do i=1,10    !loop over each parameter for a single planet
      ii=8+i+nplanet*(np-1)
      sol_wl(ii)=sol(solrange(ii,1))
      solerr_wl(ii,:)=solerr(solrange(ii,1),:)
      solrange_wl(ii,:)=ii
   enddo
enddo

nwv_wl=1 !there is only one band pass for whitelight
!Fit model to the whitelight observations
call fittransitmodel8(npars_wl,nplanet,sol_wl,solerr_wl,solrange_wl,    &
 nwv_wl,nobs,time_wl,flux_wl,ferr_wl,exptime_wl,ntt,tobs,omc)

!use whitelight solution to update global model
do i=1,8 !loop over global parameters
   do j=solrange(i,1),solrange(i,2)
      sol(j)=sol_wl(i)
   enddo
enddo
do np=1,nplanet !loop over each planet in transit model
   do i=1,10    !loop over each parameter for a single planet
      ii=8+i+nplanet*(np-1)
      do j=solrange(ii,1),solrange(ii,2)
         sol(j)=sol_wl(ii)
      enddo
   enddo
enddo

!The following part is for testing the whitelight solution
!make a transit-model to compare to the data
ntest=1 !ntest=0 - run section, otherwise, skip
if(ntest.eq.0)then
   allocate(sptmodel(nwv_wl,nobs)) !array to hold the spectral transit model
   nwvc=0 !calculate model for all bandpasses
   call sptransitmodel(nplanet,npars_wl,sol_wl,solrange_wl,nwv_wl,nobs, &
    time_wl,exptime_wl,ntt,tobs,omc,sptmodel,nwvc)

   !output to check average is good.
   do i=1,nobs
      write(6,500) time_wl(1,i),flux_wl(1,i),ferr_wl(1,i),sptmodel(1,i)
   enddo
   500 format(10000(1PE17.10,1X))
   write(0,*) "Ready..."
   read(5,*)
endif

return
end subroutine fitwhitelighttransit
