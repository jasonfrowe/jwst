subroutine modelline(npt,line,ntrace,sol,isol)
use precision
use fittingmod
implicit none
integer, target :: ntrace,nfit
integer :: npt,info,lwa,nfitin,j,i
integer, dimension(:), target :: isol
integer, allocatable, dimension(:) :: iwa
real(double) :: tol
real(double), dimension(:), target :: line,sol
real(double), allocatable, dimension(:) :: solin
real(double), allocatable, dimension(:) :: wa,fvec
external fcn

!update pointers for fittingmod parameters to be passed to fvec via mod
nfit2 => nfit     !number of parameters in model
ntrace2 => ntrace !number of traces to model
line2 => line     !the data to model
isol2 => isol     !which parameters to fit
sol2 => sol       !the fitted solution

nfit=1+9*ntrace !size of sol
tol=1.0d-8 !fitting tolerence
info=0     !keep track of errors from minpack
lwa=npt*nfit+5*npt*nfit  !work-space for fitter
allocate(wa(lwa))
allocate(fvec(npt)) !contains chi-sq vectors
allocate(iwa(nfit)) !work space for fitter
!count number of parameters to fit
allocate(solin(nfit)) !allocate space for solin
nfitin=0
do i=1,nfit  !check which variables we are fitting
   if(isol(i).ne.0)then
      nfitin=nfitin+1
      solin(nfitin)=sol(i)
   endif
enddo

!write(0,*) "Begin Fitting.."
call lmdif1(fcn,npt,nfitin,solin,fvec,tol,info,iwa,wa,lwa)
if(info.gt.3)then
   write(0,*) "info: ",info
endif

!move solin into sol
j=0
do i=1,nfit
   if(isol(i).ne.0)then
      j=j+1
      sol(i)=solin(j)
   endif
enddo

!write(0,*) "sol(1): ",sol(1)

return
end subroutine modelline

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(m,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: m,n,iflag,i,j
real(double), dimension(n) :: x
real(double), dimension(m) :: fvec
real(double), allocatable, dimension(:) :: sol

interface
   subroutine psfmodel1d(npt,model,ntrace,sol)
      use precision
      implicit none
      integer, intent(in) :: npt,ntrace
      real(double), dimension(:), intent(in) :: sol
      real(double), dimension(:), intent(inout) :: model
   end subroutine psfmodel1d
end interface

!m==npt
!x==solin
allocate(sol(nfit2))
sol=sol2 !copy over input
!move solin into sol
j=0
do i=1,nfit2
   if(isol2(i).ne.0)then
      j=j+1
      sol(i)=x(j)
   endif
enddo

!write(0,*) (sol(i),i=1,nfit2)

call psfmodel1d(m,fvec,ntrace2,sol)
!read(5,*)

do i=1,m
   fvec(i)=(fvec(i)-line2(i))/1.0d0
enddo

return
end
