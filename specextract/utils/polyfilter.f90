subroutine polyfilter(nlines,nTrace,apfluxl,apfluxu,nbin)
!Jason Rowe - jasonfrowe@gmail.com
use precision
implicit none
integer :: nlines,nTrace,nbin
real(double), dimension(:,:) :: apfluxl,apfluxu
!routine specific variables
integer, allocatable, dimension(:) :: ia
integer :: i,j,k,n1,n2,nbd2,nbinp1,ndat,nfit
real(double) :: f,chisq
real(double), allocatable, dimension(:) :: x,y,sig,filt,a
real(double), allocatable, dimension(:,:) :: covar

!parameters
nfit=3
!arrays for polynomial fitter
allocate(a(nfit),ia(nfit),covar(nfit,nfit))
ia=1

!constants
nbd2=nbin/2
nbinp1=nbin+1

!allocate arrays for fitting
allocate(x(nbinp1),y(nbinp1),sig(nbinp1))
sig=0.0001 !used for weighted fittings
!allocate array to hold continuum estimate
allocate(filt(nlines))

k=1 !Trace to work on
!filter lower trace
do i=1,nlines
   if(apfluxl(i,k).gt.0.0d0)then
      n1=max(1,i-nbd2)
      n2=min(nlines,i+nbd2)
      ndat=0
      do j=n1,n2
         if(apfluxl(j,k).gt.0.0d0)then
            ndat=ndat+1
            x(ndat)=dble(j-i)
            y(ndat)=apfluxl(j,k)
         endif
      enddo
      if(ndat.gt.nfit+1)then !make sure we have enough data for a fit
         call lfit(x,y,sig,ndat,a,ia,nfit,covar,nfit,chisq)
         f=a(1)
         filt(i)=f
      elseif(ndat.eq.0)then !no data.. no calculations needed
         filt(i)=0.0d0
      else !use average if we cannot fit data
         filt(i)=Sum(y(1:ndat))/dble(ndat)
      endif
   else
      filt(i)=0.0d0 !avoid corrections when flux=0
   endif
enddo
apfluxl(1:nlines,k)=apfluxl(1:nlines,k)-filt(1:nlines)

!filter upper trace
do i=1,nlines
   if(apfluxu(i,k).gt.0.0d0)then
      n1=max(1,i-nbd2)
      n2=min(nlines,i+nbd2)
      ndat=0
      do j=n1,n2
         if(apfluxu(j,k).gt.0.0d0)then
            ndat=ndat+1
            x(ndat)=dble(j-i)
            y(ndat)=apfluxu(j,k)
         endif
      enddo
      if(ndat.gt.nfit+1)then !make sure we have enough data for a fit
         call lfit(x,y,sig,ndat,a,ia,nfit,covar,nfit,chisq)
         f=a(1)
         filt(i)=f
      elseif(ndat.eq.0)then !no data.. no calculations needed
         filt(i)=0.0d0
      else !use average if we cannot fit data
         filt(i)=Sum(y(1:ndat))/dble(ndat)
      endif
   else
      filt(i)=0.0d0 !avoid corrections when flux=0
   endif
enddo
apfluxu(1:nlines,k)=apfluxu(1:nlines,k)-filt(1:nlines)

return
end subroutine polyfilter
