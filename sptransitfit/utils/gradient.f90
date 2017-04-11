subroutine gradient(nwv,nobs,nplanet,npars,sol,solerr,solrange,time,    &
 flux,ferr,exptime,ntt,tobs,omc,f,g)
use precision
implicit none
!import vars
integer :: nwv,nobs,nplanet,npars
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double) :: f
real(double), dimension(:) :: sol,g
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc
!local vars
integer :: i,j,inpar
integer, allocatable, dimension(:) :: ifpar
real(double) :: small,ftest1
real(double), allocatable, dimension(:) :: soltest

interface
   function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,     &
    flux,ferr,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: time,flux,ferr,exptime,  &
       tobs,omc
      real(double) :: loglikelihood
   end function loglikelihood
end interface

allocate(soltest(npars),ifpar(npars)) !array to hold preturbed model solution

small=1.0d-8
j=0 !counter for number of fitted parameters
do i=1,npars
   if(solerr(i,1).ne.0.0)then
      j=j+1 !increase counter to track number of fitted parameters
      ifpar(j)=i !precompute with parameters have gradient calc. to paral.
   endif
enddo
inpar=j

!now we calculate gradient for each fitted parameter.  This task is
!highly parallelizable.
!$OMP PARALLEL DO
do j=1,inpar
   i=ifpar(j)
   soltest=sol !make copy of original solution
   soltest(i)=soltest(i)+small !add small perturbation
   ftest1=-loglikelihood(nwv,nobs,nplanet,npars,soltest,solrange,       &
    time,flux,ferr,exptime,ntt,tobs,omc)
   g(j)=(ftest1-f)/small
enddo
!$OMP END PARALLEL DO

return
end subroutine gradient
