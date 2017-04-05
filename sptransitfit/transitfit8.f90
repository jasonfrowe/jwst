program transitfit8
!(c) Jason Rowe 2017 jason.rowe@ubishops.ca
use precision
implicit none
integer iargc,nobsmax,nwvmax,nunit,nobs,nwv,npars,nparsmax,i,j,         &
 nplanetmax,nplanet
real(double), allocatable, dimension(:) :: time,sol
real(double), allocatable, dimension(:,:) :: flux,solerr
character(80) :: obsfile,parsfile

interface
   subroutine getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,   &
    solerr)
      use precision
      implicit none
      integer :: nparsmax,npars,nunit,nplanet,nplanetmax
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr
   end subroutine getfitpars
end interface

if(iargc().lt.1)then
   write(6,*) "Usage: transitfit5 <specfile> <modelpars>"
   write(6,*) "  <specfile>  : file containing extracted time-series spectra"
   write(6,*) "  <modelpars> : file containing model parameters for fitting"
   stop
endif

!get filename for input data file from first commandline arguement
call getarg(1,obsfile)
!get filename for input model parameters
call getarg(2,parsfile)

!opens file and assigns unit number.  source in readdata.f90
call openfile(nunit,obsfile)
!scan through the observations and read in the number of bandpasses
! and observations.  Source in readdata.f90
call getnobsnwv(nunit,nobsmax,nwvmax)
!write(0,*) "Maximum number of bandpasses:   ",nwvmax
!write(0,*) "Maximum number of observations: ",nobsmax

!now we can allocate arrays to contain observation times and fluxes
allocate(time(nobsmax),flux(nwvmax,nobsmax))

!read in data
call readdata(nunit,nobsmax,nwvmax,nobs,nwv,time,flux)
write(0,*) "Number of bandpasses:   ",nwv
write(0,*) "Number of observations: ",nobs

!close data file
close(nunit)

!Now we move on to reading in model parameters.

!opens file and assigns unit number.  source in readdata.f90
call openfile(nunit,parsfile)

!calculate the number of parameters we are fitting
call getnumfitpars(nunit,nparsmax,nplanetmax)
!write(0,*) "Maximum number of model parameters: ",nparsmax
!write(0,*) "Maximum number of planets in model: ",nplanetmax

!allocate variables for model parameters
allocate(sol(nparsmax),solerr(nparsmax,4))

!read in model parameters
call getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,solerr)
write(0,*) "Number of model parameters read: ",npars
write(0,*) "Number of planets in model: ",nplanet
!do i=1,npars
!   write(6,500) sol(i),(solerr(i,j),j=1,4)
!   500 format(40(1X,1PE17.10))
!enddo

!close parameter file
close(nunit)

end program transitfit8
