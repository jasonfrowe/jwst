program sptransitplot8
use precision
implicit none
integer :: iargc,nunit,nobsmax,nwvmax,nobs,nwv,nparsmax,nplanetmax,     &
 npars,nplanet,i,nwvc
integer, allocatable, dimension(:) :: ntt
integer, allocatable, dimension(:,:) :: solrange
real(double) :: chisq,norm1,resfrac
real(double),allocatable, dimension(:) :: sol,rdr,norm
real(double), allocatable, dimension(:,:) :: time,flux,ferr,exptime,    &
 solerr,tobs,omc,sptmodel,res
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
   subroutine plotimg(nwv,nobs,res,resfrac)
      use precision
      implicit none
      integer :: nwv,nobs
      real(double) :: resfrac
      real(double), dimension(:,:) :: res
   end subroutine plotimg
   subroutine plotrdr(nwv,rdr)
      use precision
      implicit none
      integer :: nwv
      real(double), dimension(:) :: rdr
   end subroutine plotrdr
   subroutine plotres(nwv,nobs,res)
      use precision
      implicit none
      integer :: nwv,nobs
      real(double),dimension(:,:) :: res
   end subroutine plotres
end interface

if(iargc().lt.1)then
   write(6,*) "Usage: sptransitplot8 <specfile> <modelpars> <ttfiles>"
   write(6,*) "  <specfile>  : file containing extracted time-series spectra"
   write(6,*) "  <modelpars> : file containing model parameters for fitting"
   write(6,*) "  <ttfiles>   : file(s) containing TTVs. One for each planet"
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
!write(0,*) "Maximum number of model parameters: ",nparsmax
!write(0,*) "Maximum number of planets in model: ",nplanetmax

!allocate variables for model parameters
!sol contains the nominal values of the model
!solerr(,1) dicates whether a parameter is fit. 0=no.
!           Also controls Gibbs for MCMC
!solerr(,2) is mode for prior
!solerr(,3) is +1sig for prior
!solerr(,4) is -1sig for prior (value should be negative)
!solrange() gives inclusive range of indices for specific parameter
!           1:RHO,2:NL1,3:NL2,4:NL3,5:NL4,6:DIL,7:VOF,8:ZPT,9:EP,
!           10:PE,11:BB,12:RD,13:EC,14:ES,15:KR,16:TE,17:EL,18:AL
allocate(sol(nparsmax),solerr(nparsmax,4),solrange(8+10*nplanetmax,2))

!read in model parameters
call getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,solerr,     &
 solrange)
write(0,*) "Number of model parameters read: ",npars
write(0,*) "Number of planets in model: ",nplanet

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

!close parameter file
close(nunit)

!open PGPLOT device
call pgopen('?')
call pgsubp(2,2)
call pgpage()
call PGPAP (12.0 ,0.5) !use a rectangle 8" long
call pgslw(3) !thicker lines
call pgsch(2.0) !bigger text

allocate(sptmodel(nwv,nobs)) !array to hold the spectral transit model
nwvc=0
call sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,   &
   ntt,tobs,omc,sptmodel,nwvc)

!print out reduced chi-sq
chisq=Sum((sptmodel-flux)**2.0d0/(ferr*ferr))
write(0,*) "reduced chi-sq: ",chisq/dble(nwv*nobs)

!get normalization
if (solrange(8,1).eq.solrange(8,2)) then
   allocate(norm(1))
   norm1=sol(solrange(8,1))
else
   allocate(norm(nwv))
   norm(1:nwv) = sol(solrange(8,1):solrange(8,2))
endif

!plot image of residuals
!calculate residuals
allocate(res(nwv,nobs))
res=(flux-sptmodel) !residuals

!normalize residuals
if(size(norm).gt.1)then
   do i=1,nwv
      res(i,:)=res(i,:)/norm(i)
   enddo
else
   res=res/norm1
endif

resfrac=0.1 !range of data to display, must be between 0 and 1.0, 1.0 is to use all data
call pgpanl(2,2)
call plotimg(nwv,nobs,res,resfrac)

!alternate residual plot
call pgpanl(1,2)
call plotres(nwv,nobs,res)

!plot r/R*
!get time-series of RDR
allocate(rdr(nwv))
rdr=sol(solrange(12,1):solrange(12,2))
call pgpanl(1,1)
call plotrdr(solrange(12,2)-solrange(12,1)+1,rdr)

call pgclos()

end program sptransitplot8
