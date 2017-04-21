subroutine gradient(nwv,nobs,nplanet,npars,sol,solerr,solrange,time,    &
 flux,ferr,exptime,ntt,tobs,omc,f,g,sptmodel)
use precision
implicit none
!import vars
integer :: nwv,nobs,nplanet,npars
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double) :: f
real(double), dimension(:) :: sol,g
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc, &
 sptmodel
!local vars
integer :: i,ii,j,k,kk,l,inpar,np,nwvc
integer, allocatable, dimension(:) :: ifpar,inwv
real(double) :: small,ftest1,h
real(double), allocatable, dimension(:) :: soltest
real(double), allocatable, dimension(:,:) :: sptmodel1

interface
   function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,     &
    flux,ferr,exptime,ntt,tobs,omc,sptmodel,nwvc)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars,nwvc
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: time,flux,ferr,exptime,  &
       tobs,omc,sptmodel
      real(double) :: loglikelihood
   end function loglikelihood
end interface

allocate(soltest(npars),ifpar(npars),inwv(npars)) !array to hold preturbed model solution

small=1.5d-8

j=0 !number of fitted parameters
l=0 !counter as we walk though model parameters
do i=1,8
   kk=0
   do k=solrange(i,1),solrange(i,2)
      l=l+1
      kk=kk+1 !increase counter to track bandpasss
      if(solerr(l,1).ne.0.0)then
         j=j+1 !increase counter to track number of fitted parameters
         ifpar(j)=l !precompute with parameters have gradient calc. to paral.
         if((solrange(i,2)-solrange(i,1)).eq.0)then
            inwv(j)=0  !parameter applied to all bandpasses
         else
            inwv(j)=kk  !parameter only applies to a single bandpass
         endif
      endif
   enddo
enddo

!planet parameters
do np=1,nplanet !loop over each planet in transit model
   do i=1,10    !loop over each parameter for a single planet
      ii=8+i+nplanet*(np-1)
      kk=0 !keep track of bandpass
      do k=solrange(ii,1),solrange(ii,2)
         l=l+1
         kk=kk+1 !increase counter to track bandpasss
         if(solerr(l,1).ne.0.0)then
            j=j+1
            ifpar(j)=l
            if((solrange(ii,2)-solrange(ii,1)).eq.0)then
               inwv(j)=0  !parameter applied to all bandpasses
            else
               inwv(j)=kk  !parameter only applies to a single bandpass
            endif
         endif
      enddo
   enddo
enddo
inpar=j

!used to make new likelihood model
allocate(sptmodel1(nwv,nobs))
sptmodel1=0.0d0

!now we calculate gradient for each fitted parameter.  This task is
!highly parallelizable.
!$OMP PARALLEL DO PRIVATE(i,soltest,ftest1,nwvc,sptmodel1,h)
do j=1,inpar
   i=ifpar(j)
   soltest=sol !make copy of original solution
   h=small*dabs(soltest(i))
   if((h.eq.0.0d0).or.(soltest(i).eq.soltest(i)+h)) h=small
   soltest(i)=soltest(i)+h !add small perturbation
   nwvc=inwv(j)
   sptmodel1=sptmodel
   ftest1=-loglikelihood(nwv,nobs,nplanet,npars,soltest,solrange,       &
    time,flux,ferr,exptime,ntt,tobs,omc,sptmodel1,nwvc)
   g(j)=(ftest1-f)/h
   !write(0,501) i,sol(i),soltest(i),h,g(j)
   !read(5,*)
   !if(j.eq.62) write(0,500) j,i,sol(i),soltest(i)
enddo
!$OMP END PARALLEL DO
500 format(2(I3,1X),3(1PE17.10,1X))
501 format(I2,1X,4(1PE17.10,1X))

return
end subroutine gradient
