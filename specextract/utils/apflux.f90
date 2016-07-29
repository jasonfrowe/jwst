subroutine apflux(naxes,Image,bpix,nTrace,dTrace,nap,nsky,flux)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nTrace,i,nap,nsky,k,nm,np,nsm,nsp
integer, dimension(2) :: naxes
real(double) :: bpix,sky
real(double), dimension(:,:) :: Image,dTrace
real(double), dimension(:) :: flux
!plotting
!real :: pmin
!real, allocatable, dimension(:) :: px,py,bb

!nap=20 !aperture for flux sums
!nsky=50 !aperture for sky estimate

k=1 !working on first trace for now.
do i=1,naxes(1)
   !range of pixels for flux
   nm=max(1,int(dTrace(i,k))-nap)
   np=min(naxes(2),int(dTrace(i,k))+nap)
   !range of pixels for sky
   nsp=min(naxes(2),np+nsky+1)
   nsm=max(1,nm-nsky-1)
   sky=0.0d0
   if(nsp.gt.np+1) sky=sky+sum(Image(i,np+1:nsp))
   if(nsm.lt.nm-1) sky=sky+sum(Image(i,nsm:nm-1))
   flux(i)=Sum(Image(i,nm:np))-sky
enddo

!open(unit=11,file='apflux.dat')
!do i=1,naxes(1)
!   write(11,*) i,flux(i)
!enddo
!close(11)



return
end subroutine apflux
