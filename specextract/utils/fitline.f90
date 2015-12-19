subroutine fitline(nlines,nTrace,apfluxl,apfluxu,nfit)
!Jason Rowe jasonfrowe@gmail.com
use precision
implicit none
integer :: nfit,nlines,nTrace
real(double), dimension(:,:) :: apfluxl,apfluxu

integer :: i,j,k,ndat
integer, allocatable, dimension(:) :: ia
real(double) :: chisq,f
real(double), allocatable, dimension(:) :: a,xf,yf,sig
real(double), allocatable, dimension(:,:) :: covar

k=1
allocate(xf(nlines),yf(nlines),sig(nlines))
ndat=0
do i=1,nlines
   if(apfluxl(i,k).gt.0.0)then
      ndat=ndat+1
      xf(ndat)=dble(i)
      yf(ndat)=apfluxl(i,k)
   endif
enddo
sig=0.0001 !ideally this should be set to the standard deviation
allocate(a(nfit),ia(nfit),covar(nfit,nfit))
ia=1 !fit all variables
call lfit(xf,yf,sig,ndat,a,ia,nfit,covar,nfit,chisq)
!write(0,*) "yo:",a
do i=1,nlines
   f=a(1)
   do j=2,nfit
      f=f+a(j)*dble(i)**dble(j-1)
   enddo
   apfluxl(i,k)=apfluxl(i,k)-f
enddo
ndat=0
do i=1,nlines
   if(apfluxu(i,k).gt.0.0)then
      ndat=ndat+1
      xf(ndat)=dble(i)
      yf(ndat)=apfluxu(i,k)
   endif
enddo
call lfit(xf,yf,sig,ndat,a,ia,nfit,covar,nfit,chisq)
!write(0,*) "yo:",a
do i=1,nlines
   f=a(1)
   do j=2,nfit
      f=f+a(j)*dble(i)**dble(j-1)
   enddo
   apfluxu(i,k)=apfluxu(i,k)-f
enddo
deallocate(xf,yf,sig,a,ia,covar)

return
end subroutine fitline
