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
integer :: i,j,nld1,nld(4)
real(double) Pi,tPi,ll1,ll2,ll3,ld(4),prior
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
! *** valid limb-darkening *** !
nld1=0
nld=0
do i=1,4
   nld(i)=solrange(i+1,2)-solrange(i+1,1)+1
   nld1=max(nld(i),nld1) !how many LD pars vary (max)
enddo
do i=1,nld1
   do j=1,4
      if(nld(j).eq.nld1)then
         ld(j)=sol(solrange(j+1,1)+i-1)
      else
         ld(j)=sol(solrange(j+1,1))
      endif
   enddo

   prior=1.0
   if((ld(1).eq.0.0).and.(ld(2).eq.0.0))then
      if((ld(3).lt.0.0).or.(ld(3).gt.1.0).or.(ld(4).lt.0.0).or.         &
       (ld(4).gt.1.0))then
         prior=9.9d30
      endif
   elseif((ld(3).eq.0.0).and.(ld(4).eq.0.0))then
   !write(0,*) "Quad..",c1,c2
      if((ld(1)+ld(2).gt.1.0).or.(ld(1).lt.0).or.                       &
       (ld(1)+2.0d0*ld(2).lt.0))then
         !write(0,*) "invalid.."
         prior=9.9d30
      endif
   endif
!   write(0,500) ld,prior
enddo
!write(0,*) "prior: ",prior
!read(5,*)

500 format(10000(1PE17.10,1X))

! *** 0 < e < 1 *** !

!apply prior
loglikelihood=loglikelihood*prior

return
end function loglikelihood
