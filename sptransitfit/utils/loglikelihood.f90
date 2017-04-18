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

!check model constraints, no point calculating a model if parameter
!is out of bounds

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

   if((ld(1).eq.0.0).and.(ld(2).eq.0.0))then
      if((ld(3).lt.0.0).or.(ld(3).gt.1.0).or.(ld(4).lt.0.0).or.         &
       (ld(4).gt.1.0))then
         loglikelihood=9.9d30
         return
      endif
   elseif((ld(3).eq.0.0).and.(ld(4).eq.0.0))then
   !write(0,*) "Quad.."
      if((ld(1)+ld(2).gt.1.0).or.(ld(1).lt.0).or.                       &
       (ld(1)+2.0d0*ld(2).lt.0))then
         !write(0,*) "invalid.."
         loglikelihood=9.9d30
         return
      endif
   endif
enddo


500 format(10000(1PE17.10,1X))

! *** 0 < e < 1 *** !

!make a transit-model to compare to the data
allocate(sptmodel1(nwv,nobs))
call sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,   &
   ntt,tobs,omc,sptmodel1,nwvc)
if(nwvc.eq.0)then
   sptmodel=sptmodel1
else
   sptmodel(nwvc,:)=sptmodel1(nwvc,:)
endif


!begin calculating log loglikelihood
!this assumes all data is valid.
ll1=npars*nwv*log(tPi) !number of observations

allocate(ferr2(nwv,nobs))
ferr2=ferr*ferr
ll2=Sum(log(ferr2)) !vectorized (seems to be faster)

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
