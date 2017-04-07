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
real(double), allocatable, dimension(:,:) :: bt

!deal with case of no out-of-transit data.

interface
   subroutine getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,   &
    omc,bt)
      use precision
      implicit none
      integer :: nwv,nplanet,npars,nobs
      integer, dimension(:) :: ntt
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: bt,solrange,time,tobs,omc
   end subroutine getb
end interface

!calculate if data is in transit or out of transit.
!this means calculating b and using rprs.
allocate(bt(nwv,nobs))
write(0,*) "Calling getb"
call getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,omc,bt)

!check if we have one zpt per bandpass, or just one.
if(solrange(8,2)-solrange(8,1).gt.0)then !8 corresponds to zero-point

   !gather all observations that are out of transit and calculate median

else !else calculate for each bandpass

   !on a wavelength by wavelength basis get out of tranist data and
   !calculate median

endif

return
end subroutine estzpt
