program transitmcmc8
use precision
implicit none
integer :: nunit,nobsmax,nwvmax,nobs,nwv,nparsmax,nplanetmax,npars,nplanet,i, &
 seed,nwvc
integer, allocatable, dimension(:) :: ntt
integer, allocatable, dimension(:,:) :: solrange
real(double) :: logl1
real(double), allocatable, dimension(:) :: sol
real(double), allocatable, dimension(:,:) :: time,flux,ferr,exptime,solerr, &
 tobs,omc,sptmodel
character(80) :: obsfile,parsfile,ttfile

interface
   subroutine getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,   &
    solerr,solrange)
      use precision
      implicit none
      integer :: nparsmax,npars,nunit,nplanet,nplanetmax
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr
   end subroutine getfitpars
   subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,  &
      exptime,ntt,tobs,omc,sptmodel,nwvc)
      use precision
      implicit none
      integer :: nplanet,npars,nwv,nobs,nwvc
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: sptmodel,tobs,omc,time,exptime
   end subroutine sptransitmodel
   subroutine mhgmcmc(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,ntt,tobs,&
    omc,seed)
      use precision
      implicit none
      integer :: nplanet,npars,nwv,nobs,seed
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: tobs,omc,time,exptime
   end subroutine mhgmcmc
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
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: transitmcmc8 <specfile> <modelpars> <ttfiles> <niter>"
   write(0,*) "  <specfile>  : file containing extracted time-series spectra"
   write(0,*) "  <modelpars> : file containing model parameters for fitting"
   write(0,*) "  <ttfiles>   : file(s) containing TTVs. One for each planet"
   write(0,*) "    <niter>   : number of MCMC iterations (Chain length)"
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

!now we can allocate arrays to contain observation times and fluxes
allocate(time(nwvmax,nobsmax),flux(nwvmax,nobsmax),ferr(nwvmax,nobsmax),&
 exptime(nwvmax,nobsmax))

!read in data
call readdata(nunit,nobsmax,nwvmax,nobs,nwv,time,flux,ferr,exptime)
write(0,*) "Number of bandpasses:   ",nwv
write(0,*) "Number of observations: ",nobs

!close data file
close(nunit)

!Now we move on to reading in model parameters.

!opens file and assigns unit number.  source in readdata.f90
call openfile(nunit,parsfile)

!calculate the number of parameters we are fitting
call getnumfitpars(nunit,nparsmax,nplanetmax)

!allocate arrays to hold model parameters
allocate(sol(nparsmax),solerr(nparsmax,4),solrange(8+10*nplanetmax,2))

!read in model parameters
call getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,solerr,     &
 solrange)
write(0,*) "Number of model parameters read: ",npars
write(0,*) "Number of planets in model: ",nplanet

!close parameter file
close(nunit)

!read in transit timing variations
!ntt contains number of TT measurements, tobs is the observer location
!of the transit and OMC is the observed minus calculated measurement
allocate(ntt(nplanet),tobs(nplanet,nobs),omc(nplanet,nobs))
do i=1,nplanet
   if(iargc().ge.2+i)then
      call getarg(2+i,ttfile)
      if(ttfile.eq.'null')then
         ntt(i)=0
      else
         nunit=10
         open(unit=nunit,file=ttfile,status='old',err=905)
         goto 906
          905 write(0,*) "Cannot open ", ttfile
          stop
         906 continue
         call readttfile(nunit,nplanet,nobs,i,ntt,tobs,omc)
         close(nunit)
      endif
   else
      ntt(i)=0
   endif
enddo

!initialize random number generator
call initran(seed)

!log-likelihood is calculated for all bandpasses
nwvc=0
!allocate sptmodel to hold global model solution
allocate(sptmodel(nwv,nobs))
!calculate likelihood for nominal model
logl1=-loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,flux,   &
 ferr,exptime,ntt,tobs,omc,sptmodel,nwvc)

write(0,*) "logl1: ",logl1

!call M-H-G to get buffer.
call mhgmcmc(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,ntt,tobs, &
 omc,seed)




end program transitmcmc8

