subroutine fittransitmodel8v2(npars,nplanet,sol,solerr,solrange,nwv,    &
 nobs,time,flux,ferr,exptime,ntt,tobs,omc)
use precision
use fittermod
implicit none
!import vars
integer, target :: npars,nwv,nobs,nplanet
integer, dimension(:), target :: ntt
integer, dimension(:,:), target :: solrange
real(double), dimension(:), target :: sol
real(double), dimension(:,:), target :: solerr,time,flux,ferr,exptime,  &
 tobs,omc
!local vars
integer :: n,i
integer, allocatable, dimension(:) :: ia
real(double) :: mlogl,alambda,mloglgoal
real(double), allocatable, dimension(:) :: atry,beta,da
real(double), allocatable, dimension(:,:) :: alpha,covar
!export fit
integer :: nunit,filestatus
character(80) :: newfitfile


interface
   subroutine EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,   &
    flux,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: npars,nwv,nobs,nplanet
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc
   end subroutine EstZpt
   subroutine lmmin(ma,a,ia,mlogl,alambda,atry,beta,da,alpha,covar,mloglgoal)
      use precision
      implicit none
      integer :: ma
      integer, dimension(:) :: ia
      real(double) :: mlogl,alambda,mloglgoal
      real(double), dimension(:) :: a,atry,beta,da
      real(double), dimension(:,:) :: alpha,covar
   end subroutine lmmin
   subroutine exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
      use precision
      implicit none
      integer :: nunit,npars,nplanet
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr
   end subroutine exportfitpars
end interface

!Get estimate for zero points.
call EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,flux,       &
 exptime,ntt,tobs,omc)

!pick off variables that are being fitted.
allocate(ia(npars))
ia=0 !default is to fit nothing
do i=1,npars
   if(solerr(i,1).ne.0.0)then
      ia(i)=1 !flag variable to be fit
   endif
enddo

!allocate work arrays
allocate(atry(npars),beta(npars),da(npars),alpha(npars,npars),          &
 covar(npars,npars))

!update pointers
nwv2 => nwv
nobs2 => nobs
nplanet2 => nplanet
npars2 => npars
!sol2 => sol
solerr2 => solerr
solrange2 => solrange
time2 => time
flux2 => flux
ferr2 => ferr
exptime2 => exptime
ntt2 => ntt
omc2 => omc

!ideal mlogl
!mloglgoal=-(nobs*nwv*1.837877+Sum(log(ferr*ferr))+nobs*nwv)
mloglgoal=nobs*nwv
!write(0,*) "ideal mlogl: ",mloglgoal

alambda=-1.0
!write(0,*) "First call to lmmin"
call lmmin(npars,sol,ia,mlogl,alambda,atry,beta,da,alpha,covar,mloglgoal)
!write(0,*) "done.."
do while(alambda.lt.1.0d8) !i=1,1000
   alambda=max(alambda,1.0d-10)
   call lmmin(npars,sol,ia,mlogl,alambda,atry,beta,da,alpha,covar,mloglgoal)
!   write(0,*) "RD: ",sol(12),mlogl/mloglgoal
!   read(5,*)

!   newfitfile="newfit.dat"
!   nunit=10
!   open(unit=nunit,file=newfitfile,iostat=filestatus)
!   if(filestatus>0)then !trap missing file errors
!      write(0,*) "Cannot open ",newfitfile
!      stop
!   endif
!   call exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
!   close(nunit)

enddo

return
end subroutine fittransitmodel8v2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine lmmin(ma,a,ia,mlogl,alambda,atry,beta,da,alpha,covar,mloglgoal)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!heavily influnced by mrqmin from Numerical Recipes
!set alambda < 0 for initialization
!ma - number of parameters in model (fitted and non-fitted)
!a - array of model parameters (a(ma))
!ia - int array to flag which variables are fit. zero=fixed, otherwise fit
!alambda - levenberg-marquardt fudge factor
use precision
implicit none
!import vars
integer :: ma
integer, dimension(:) :: ia
real(double) :: mlogl,alambda,mloglgoal
real(double), dimension(:) :: a,atry,beta,da
real(double), dimension(:,:) :: alpha,covar
!local vars
integer :: j,k,l,info,ipiv
integer, save :: mfit
real(double), save :: omlogl
real(double), allocatable, dimension(:) :: dmloglda

interface
   subroutine lmcof(ma,mfit,a,ia,alpha,beta,mlogl,dmloglda)
      use precision
      implicit none
      integer :: ma,mfit
      integer, dimension(:) :: ia
      real(double) :: mlogl
      real(double), dimension(:) :: a,beta,dmloglda
      real(double), dimension(:,:) :: alpha
   end subroutine lmcof
   subroutine plotimg(nwv,nobs,res)
      use precision
      implicit none
      integer :: nwv,nobs
      real(double), dimension(:,:) :: res
   end subroutine plotimg
end interface

allocate(dmloglda(ma))

if(alambda.lt.0.0d0)then
   mfit=0
   do j=1,ma
      if(ia(j).ne.0) mfit=mfit+1 !count number of fitted variables
   enddo
   alambda=0.001
   call lmcof(ma,mfit,a,ia,alpha,beta,mlogl,dmloglda)
   omlogl=mlogl
   atry=a
endif

covar=alpha
da=beta
do j=1,mfit
   covar(j,j)=alpha(j,j)*(1.0d0+alambda)
enddo

call gaussj(covar,mfit,ma,da,1,1) !solve matrix

!call pgopen('/xwindow')
!call PGPAP (8.0 ,1.0) !use a square 8" across
!call pgslw(3) !thicker lines
!call plotimg(mfit,mfit,covar)
!call pgclos()

if(alambda.eq.0.d0)then
   call covsrt(covar,ma,ma,ia,mfit)
   call covsrt(alpha,ma,ma,ia,mfit)
   return
endif

!apply new solution
j=0
do l=1,ma
   if(ia(l).ne.0)then
      j=j+1
      atry(l)=a(l)+da(j)/alambda*mlogl/(mloglgoal*1.0)
   endif
enddo

call lmcof(ma,mfit,atry,ia,covar,da,mlogl,dmloglda)

!write(0,501) alambda,omlogl,mlogl,omlogl-mlogl
501 format(4(1PE17.10,1X))

if(mlogl.lt.omlogl)then
!   write(0,*) "success.."
   alambda=0.1*alambda !if success, accept the new solution
   omlogl=mlogl
   alpha=covar
   beta=da
   a=atry
else
   alambda=10.0d0*alambda
   mlogl=omlogl
endif

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine lmcof(ma,mfit,a,ia,alpha,beta,mlogl,dmloglda)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: ma,mfit
integer, dimension(:) :: ia
real(double) :: mlogl
real(double), dimension(:) :: a,beta,dmloglda
real(double), dimension(:,:) :: alpha
!local vars
integer j,k,l,m

interface
   subroutine funcs(sol,mlogl,dmloglda)
      use precision
      use fittermod
      implicit none
      real(double) :: mlogl
      real(double), dimension(:) :: sol,dmloglda
   end subroutine funcs
end interface

!funcs only calculates dmloglda for variables that are fitted.
call funcs(a,mlogl,dmloglda) !get logL and gradient
!construct alpha and beta
do j=1,mfit
   do k=1,j
      alpha(j,k)=0.5*dmloglda(j)*dmloglda(k)
   enddo
   beta(j)=-0.5*dmloglda(j)
enddo

do j=2,mfit
   do k=1,j-1
      alpha(k,j)=alpha(j,k)
   enddo
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine funcs(sol,mlogl,dmloglda)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittermod
implicit none
!import vars
real(double) :: mlogl
real(double), dimension(:) :: sol,dmloglda
!local vars
integer :: nwvc
real(double), allocatable, dimension(:,:) :: sptmodel

interface
   function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,     &
    flux,ferr,exptime,ntt,tobs,omc,sptmodel,nwvc)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars,nwvc
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: time,flux,ferr,exptime,tobs,omc,  &
       sptmodel
      real(double) :: loglikelihood
   end function loglikelihood
   subroutine gradient(nwv,nobs,nplanet,npars,sol,solerr,solrange,time, &
    flux,ferr,exptime,ntt,tobs,omc,f,g,sptmodel)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double) :: f
      real(double), dimension(:) :: sol,g
      real(double), dimension(:,:) :: solerr,time,flux,ferr,   &
       exptime,tobs,omc,sptmodel
   end subroutine gradient
end interface

allocate(sptmodel(nwv2,nobs2))
nwvc=0
mlogl=-loglikelihood(nwv2,nobs2,nplanet2,npars2,sol,solrange2,time2,   &
 flux2,ferr2,exptime2,ntt2,tobs2,omc2,sptmodel,nwvc)
call gradient(nwv2,nobs2,nplanet2,npars2,sol,solerr2,solrange2,time2,  &
 flux2,ferr2,exptime2,ntt2,tobs2,omc2,mlogl,dmloglda,sptmodel)

return
end
