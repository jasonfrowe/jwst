subroutine modelline(npt,line,ntrace,sol,isol,ilinkpsf,nline)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
use fittingmod
implicit none
integer, target :: ntrace,nfit,ilinkpsf,nline
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
ilinkpsf2 => ilinkpsf !are the PSF models linked
nline2 => nline   !which column are we fitting

nfit=1+9*ntrace !size of sol
tol=1.0d-10 !fitting tolerence
info=0     !keep track of errors from minpack
lwa=npt*nfit+5*npt*nfit  !work-space for fitter
allocate(wa(lwa))
allocate(fvec(npt)) !contains chi-sq vectors
allocate(iwa(nfit)) !work space for fitter
!count number of parameters to fit
allocate(solin(nfit)) !allocate space for solin

!are we linking PSF models across orders? If yet, then set isol as needed
if(ilinkpsf.gt.0)then
   do i=2,ntrace
      do j=1,9 !number of parameters
         if((j.ne.7).and.(j.ne.8))then !position and amplitude is still fit
            isol(j+1+9*(i-1))=0
!            write(0,*) "fixed: ",j+1+9*(i-1)
         endif
      enddo
   enddo
endif
!write(0,*) "isol: "
!do i=1,size(isol)
!   write(0,*)i,isol(i)
!enddo
!read(5,*)

nfitin=0
do i=1,nfit  !check which variables we are fitting
   if(isol(i).ne.0)then
      nfitin=nfitin+1
      solin(nfitin)=sol(i)
   endif
enddo

!write(0,*) "Begin Fitting..",nfitin
call lmdif1(fcn,npt,nfitin,solin,fvec,tol,info,iwa,wa,lwa)
if(info.gt.3)then
   write(0,*) "info: ",info
endif
!write(0,*) "Done fitting.."

!move solin into sol
j=0
do i=1,nfit
   if(isol(i).ne.0)then
      j=j+1
      sol(i)=solin(j)
   endif
enddo

!if PSFs are linked then update model to reflect link
if(ilinkpsf.gt.0)then
   do i=2,ntrace
      do j=1,9 !number of parameters
         if((j.ne.7).and.(j.ne.8))then !position and amplitude is still fit
            sol(j+1+9*(i-1))=sol(j+1)
         endif
      enddo
   enddo
endif

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

!write(0,*) "fcn"

!m==npt
!x==solin
allocate(sol(nfit2))
sol=sol2 !copy over input (guess) solution
!move solin into sol
j=0
do i=1,nfit2
   if(isol2(i).ne.0)then
      j=j+1
      sol(i)=x(j)
   endif
enddo

!write(0,*) "ilink update",size(sol),size(sol2)
!if PSFs are linked then update model to reflect link
if(ilinkpsf2.gt.0)then
   do i=2,ntrace2
      do j=1,9 !number of parameters
!         write(0,*) i,j,j+1+9*(i-1)
         if((j.ne.7).and.(j.ne.8))then !position and amplitude is still fit
            sol(j+1+9*(i-1))=sol(j+1)
         endif
      enddo
   enddo
endif

!write(0,*) (sol(i),i=1,nfit2)

call psfmodel1d(m,fvec,ntrace2,sol)
!read(5,*)

do i=1,m
   fvec(i)=(fvec(i)-line2(i))/1.0d0
enddo

if(abs(sol2(3)-sol(3)).gt.1.0)then
   fvec=9.9e30
endif
if(abs(sol2(6)-sol(6)).gt.1.0)then
   fvec=9.9e30
endif
if(abs(sol2(9)-sol(9)).gt.1.0)then
   fvec=9.9e30
endif


return
end
