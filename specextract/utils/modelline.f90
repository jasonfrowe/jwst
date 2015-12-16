subroutine modelline(npt,line,ntrace,sol)
use precision
use fittingmod
implicit none
integer :: npt,ntrace,info,lwa,nfit
integer, allocatable, dimension(:) :: iwa
real(double) :: tol
real(double), dimension(:) :: line,sol
real(double), allocatable, dimension(:) :: wa,fvec
external fcn

nfit=1+9*ntrace !size of sol
tol=1.0d-8 !fitting tolerence
info=0     !keep track of errors from minpack
lwa=npt*nfit+5*npt*nfit  !work-space for fitter
allocate(wa(lwa))
allocate(fvec(npt)) !contains chi-sq vectors
allocate(iwa(nfit))

call lmdif1(fcn,npt,nfit,sol,fvec,tol,info,iwa,wa,lwa)
write(0,*) "info: ",info

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

if(n.eq.2)then
   do i=1,m
      fvec(i)=1.0d0-x(1)*(1.0d0-mu2(i))-x(2)*(1.0d0-mu2(i))**2.0d0
   enddo
else
   do i=1,m
      fvec(i)=1.0d0
      do j=1,n
         fvec(i)=fvec(i)-x(j)*(1.0d0-mu2(i)**(dble(j)/2.0d0))
      enddo
   enddo
endif

do i=1,m
   fvec(i)=(fvec(i)-dIn2(i))/1.0d0
enddo

return
end
