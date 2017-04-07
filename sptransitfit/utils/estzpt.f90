subroutine EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,flux, &
exptime,ntt,tobs,omc)
use precision
implicit none
integer :: npars,nwv,nobs,nplanet
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc,      &
   solrange
!local vars
integer i,j,k
integer, allocatable, dimension(:) :: np
integer, allocatable, dimension(:,:) :: imarktrans
real(double) :: median
real(double), allocatable, dimension(:) :: ootdata
real(double), allocatable, dimension(:,:) :: bt

!deal with case of no out-of-transit data.

interface
   subroutine getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,   &
    omc,bt,imarktrans)
      use precision
      implicit none
      integer :: nwv,nplanet,npars,nobs
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: imarktrans
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: bt,solrange,time,tobs,omc
   end subroutine getb
end interface

!calculate if data is in transit or out of transit.
!this means calculating b and using rprs.
allocate(bt(nwv,nobs),imarktrans(nwv,nobs))
call getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,omc,bt,     &
 imarktrans)

!check if we have one zpt per bandpass, or just one.
if(solrange(8,2)-solrange(8,1).eq.0)then !8 corresponds to zero-point

   !gather all observations that are out of transit and calculate median
   allocate(ootdata(nwv*nobs)) !max number of out of transit data

   k=0
   do i=1,nwv
      do j=1,nobs
         if(imarktrans(i,j).eq.1)then !only use out of transit data
            k=k+1 !increase counter
            ootdata(k)=flux(i,j)
         endif
      enddo
   enddo
   !calculate median of out of transit data as estimate of zpt
   if(k.gt.0)then
      allocate(np(k))
      call rqsort(k,ootdata,np) !sort to find median value
      median=ootdata(k/2+1)
      sol(solrange(8,1))=median  !change zpt to median
   else  !if there are no observations then leave zpt alone
      median=1.0d0
   endif

else !else calculate for each bandpass

   !on a wavelength by wavelength basis get out of tranist data and
   !calculate median

   allocate(ootdata(nobs)) !storing out of transit data
   allocate(np(nobs))  !integer array used by rqsort and for getting median

   do i=1,nwv !loop over all bandpasses
      k=0 !counter for number of observations that are out of transit
      do j=1,nobs
         if(imarktrans(i,j).eq.1)then !only use out of transit data
            k=k+1 !increase counter
            ootdata(k)=flux(i,j) !store out of transit data
         endif
         if(k.gt.0)then
            call rqsort(k,ootdata,np) !sort to find median
            median=ootdata(k/2+1)
            sol(solrange(8,1)+i-1)=median  !change zpt to median
         else !if there are no observations then leave zpt alone
            median=1.0
         endif
      enddo

   enddo

endif

return
end subroutine estzpt
