program spgenV2
!Version 2 of Spec-generator.  This version does time-series and get FITS formating correct. 
!generates GR700-SOSS spectrum with 3-order traces + PSF + response
use precision
use response
!use response
implicit none
integer, dimension(3) :: funit !number of FITS I/O
character(200), dimension(3) :: fileout !name of FITS files
!file name vars
integer :: pid,onum,vnum,gnum,spseq,anumb,enum
character(8) :: detectorname,prodtype
!random number vars
integer, dimension(3) :: now
integer :: seed
real(double) :: ran2,dumr
!Kernel vars
integer :: nrK,nKs
real(double), dimension(:,:,:), allocatable :: rKernel
!orders : traces and responces
integer :: ntracemax
real(double), allocatable, dimension(:) :: yres1,yres2,yres3
!spectral model parameters
integer :: nmodel !number of model points read in
integer :: nmodelmax !guess for size of array (will be auto increased if too small)
real(double), dimension(:), allocatable :: wmod,fmod,wmod2,fmod2 !contain wavelength, flux info
real(double), dimension(:,:), allocatable :: nll,nll2 !limb-darkening co-efficients.
!local vars
integer :: i,j !counters
integer :: noversample,nunit,filestatus,nmodeltype,iargc,iflag
real(double) :: xout, yout,rv,b
character(80) :: cline,modelfile

interface
	subroutine writefitsphdu(fileout,funit)
		use precision
     	implicit none
     	integer :: funit
     	character(200), dimension(3) :: fileout
	end subroutine writefitsphdu
	subroutine readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag)
      use precision
      implicit none
      integer, intent(inout) :: nunit,nmodelmax,nmodel,iflag
      real(double), dimension(:), intent(inout) :: wmod, fmod
   end subroutine
   subroutine readatlas(nunit,nmodelmax,nmodel,wmod,fmod,nll,iflag)
      use precision
      implicit none
      integer, intent(inout) :: nunit,nmodelmax,nmodel,iflag
      real(double), dimension(:), intent(inout) :: wmod, fmod
      real(double), dimension(:,:), intent(inout) :: nll
   end subroutine
   subroutine readKernels(nrK,nK,rKernel,noversample)
      use precision
      implicit none
      integer,intent(inout) :: nrK,nK,noversample
      real(double), dimension(:,:,:), intent(inout) :: rKernel
   end subroutine
end interface

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Command line arguments
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
if(iargc().lt.3)then
   write(0,*) "Usage: spgen <specmodel> <noversample> <planetmodel>"
   write(0,*) "   <specmodel> - Atlas-9 stellar model"
   write(0,*) " <noversample> - is new sampling for Kernel (must be > 0)"
   write(0,*) " <planetmodel> - name of planet model (A, rprs)"
   write(0,*) "           [b] - impact parameter - optional (must be > 0)"
   stop
endif

if(iargc().ge.4)then
   call getarg(4,cline)
   read(cline,*) b
   if(b.lt.0.0d0)then
      write(0,*) "b must be positive"
      stop
   endif
else
   !default impact parameter
   b=2.0!0.3589
endif

noversample=1 !now a commandline-parameter
!get oversampling from commandline
call getarg(2,cline)
read(cline,*) noversample !read in noversample
if(noversample.le.0)then
   write(0,*) "noversample must be greater than zero"
   stop
endif

!read in a model spectrum
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!read in a model spectrum
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  Model Parameters 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!image dimensions
xout=2048  !dimensions for output image.
yout=256

!parameters that control the simulation
rv=0.0 !radial velocity shift (m/s)

!parameter controling modeltype
nmodeltype=2 !1=BT-Settl, 2=Atlas-9+NL limbdarkening

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Random Number Initialization
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!Initialization of random number
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!read in Kernels
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
nrK=30 !number of Kernels to readin
nKs=64*noversample !natural size of Kernels times oversampling
allocate(rKernel(nrK,nKs,nKs))
call readKernels(nrK,nKs,rKernel,noversample)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!read in response
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ntracemax=3 !number of traces used
call readresponse() !responce is returned via pointers from response.mod
!use global variables: nres,ld,res1,res2,res3 for response function
allocate(yres1(nres),yres2(nres),yres3(nres))
call spline(ld,res1,nres,1.d30,1.d30,yres1) !set up cubic splines
call spline(ld,res2,nres,1.d30,1.d30,yres2)
call spline(ld,res3,nres,1.d30,1.d30,yres3)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!read in a model spectrum
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!file naming
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
pid = 1 !programID
onum = 1 !observation number
vnum = 1 !visit number
gnum = 1 !group visit
spseq = 1 !parallel sequence. (1=prime, 2-5=parallel)
anumb = 1 !activity number
enum = 1 !exposure number
detectorname = 'NISRAPID' !convert this 
prodtype='cal'
call getfilename(pid,onum,vnum,gnum,spseq,anumb,enum,detectorname,prodtype,fileout)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!create FITS files and insert primary HDU for each output data product. 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 1 - simulation with no convolution,  Native resolution
! 2 - simulation with convolution, native resolution
! 3 - simulation with convolution, over-sampled resolution
call writefitsphdu(fileout(1),funit(1))
call writefitsphdu(fileout(2),funit(2))
call writefitsphdu(fileout(3),funit(3))


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Read in spectral model for star
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

nmodelmax=2000000 !initital guess at the number of data points
allocate(wmod(nmodelmax),fmod(nmodelmax))
allocate(nll(nmodelmax,4)) !limb-darkening co-efficients

iflag=0 !flag traces data i/o
nmodel=0   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
   if(nmodeltype.eq.1)then
      call readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag) !read in spectrum
      nll=0.0d0 !no limb-darkening
   elseif(nmodeltype.eq.2)then
      call readatlas(nunit,nmodelmax,nmodel,wmod,fmod,nll,iflag)
   endif
   if(iflag.eq.1) then !reallocate array space (we ran out of space)
      allocate(wmod2(nmodelmax),fmod2(nmodelmax),nll2(nmodelmax,4)) !allocate temp arrays
      wmod2=wmod   !copy over the data we read
      fmod2=fmod
      nll2=nll
      deallocate(wmod,fmod,nll) !deallocate data arrays
      nmodelmax=nmodelmax*2 !lets get more memory
      write(0,*) "warning, increasing nmodelmax: ",nmodelmax
      allocate(wmod(nmodelmax),fmod(nmodelmax),nll(nmodelmax,4)) !reallocate array
      do i=1,nmodelmax/2  !copy data back into data arrays
         wmod(i)=wmod2(i)
         fmod(i)=fmod2(i)
         do j=1,4
            nll(i,j)=nll2(i,j)
         enddo
      enddo
      deallocate(wmod2,fmod2,nll2) !deallocate temp arrays
      iflag=2  !set flag that we are continuing to read in data
      cycle !repeat data read loop
   endif
   exit !successively break from data read loop.
enddo
close(nunit) !close file.
write(0,*) "Number of star model points: ",nmodel  !report number of data points read.

!close the FITS file
call closefits(funit(1))
call closefits(funit(2))
call closefits(funit(3))

end program spgenV2