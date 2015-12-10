program spgen
!generates GR700-SOSS spectrum with 3-order traces + PSF + response
use precision
use response
implicit none
integer xmax,ymax,iargc,nunit,filestatus,nmodelmax,iflag,nmodel,i,j,npx,&
   npy,nK,xout,yout,noversample,nKd2,npt,nrK,nKs,Ks,ntrace,ntracemax
integer, dimension(2) :: ybounds
real(double) :: w2p,px,py,ptrace,dxmaxp1,dymaxp1,dnossq,rv,dnpt,wl,mSum,&
   fmodres,respond,awmod
real(double), allocatable, dimension(:) :: wmod,fmod,wmod2,fmod2,       &
   fmodbin,wv,yres1,yres2,yres3
real(double), allocatable, dimension(:,:) :: pixels,Kernel,cpixels,     &
   opixels,wKernel,wpixels,wcpixels
real(double), allocatable, dimension(:,:,:) :: rKernel
character(80) :: modelfile,fileout,cline

interface
   subroutine readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag)
      use precision
      implicit none
      integer, intent(inout) :: nunit,nmodelmax,nmodel,iflag
      real(double), dimension(:), intent(inout) :: wmod, fmod
   end subroutine
end interface
interface
   subroutine writefits(nxmax,nymax,parray,fileout)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray
      character(80) :: fileout
   end subroutine
end interface
interface
   subroutine binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmodel
      real(double), intent(inout) :: rv
      real(double), dimension(:), intent(inout) :: wv,wmod,fmod
      real(double), dimension(:), intent(inout) :: fmodbin
   end subroutine
end interface
interface
   subroutine readKernels(nrK,nK,rKernel,noversample)
      use precision
      implicit none
      integer,intent(inout) :: nrK,nK,noversample
      real(double), dimension(:,:,:), intent(inout) :: rKernel
   end subroutine
end interface

if(iargc().lt.2)then
   write(0,*) "Usage: spgen specmodel noversample"
   write(0,*) "  specmodel - BT-Stell.7 stellar model"
   write(0,*) " noversample - is new sampling for Kernel (must be > 0)"
   stop
endif

xout=2048  !dimensions for output image.
yout=512

noversample=1 !now a commandline-parameter
!get oversampling from commandline
call getarg(2,cline)
read(cline,*) noversample !read in noversample
if(noversample.le.0)then
   write(0,*) "noversample must be greater than zero"
   stop
endif

if(noversample.lt.1)then
   write(0,*) "noversample must be at least 1"
   stop
endif

!read in Kernels
nrK=30 !number of Kernels to readin
nKs=64*noversample !natural size of Kernels times oversampling
allocate(rKernel(nrK,nKs,nKs))
call readKernels(nrK,nKs,rKernel,noversample)

!read in response
ntracemax=3 !number of traces used
call readresponse() !responce is returned via pointers from response.mod
!use global variables: nres,ld,res1,res2,res3 for response function
allocate(yres1(nres),yres2(nres),yres3(nres))
call spline(ld,res1,nres,1.d30,1.d30,yres1) !set up cubic splines
call spline(ld,res2,nres,1.d30,1.d30,yres2)
call spline(ld,res3,nres,1.d30,1.d30,yres3)

!read in a model spectrum
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

nmodelmax=1500000 !initital guess at the number of data points
allocate(wmod(nmodelmax),fmod(nmodelmax))
iflag=0 !flag traces data i/o
nmodel=0   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
   call readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag) !read in spectrum
   if(iflag.eq.1) then !reallocate array space (we ran out of space)
      allocate(wmod2(nmodelmax),fmod2(nmodelmax)) !allocate temp arrays
      wmod2=wmod   !copy over the data we read
      fmod2=fmod
      deallocate(wmod,fmod) !deallocate data arrays
      nmodelmax=nmodelmax*2 !lets get more memory
      write(0,*) "warning, increasing nmodelmax: ",nmodelmax
      allocate(wmod(nmodelmax),fmod(nmodelmax)) !reallocate array
      do i=1,nmodelmax/2  !copy data back into data arrays
         wmod(i)=wmod2(i)
         fmod(i)=fmod2(i)
      enddo
      deallocate(wmod2,fmod2) !deallocate temp arrays
      iflag=2  !set flag that we are continuing to read in data
      cycle !repeat data read loop
   endif
   exit !successively break from data read loop.
enddo
close(nunit) !close file.
!write(0,*) "nmodel: ",nmodel  !report number of data points read.

npt=200000 !this sets the number of spectral points when resampled
allocate(wv(npt),fmodbin(npt))
!resample model spectra on a uniform grid from 1000-30000 A
dnpt=dble(npt)
do i=1,npt
   wv(i)=1000.0+(30000.0-1000.0)/dnpt*dble(i) !make a wavelength grid
!   write(0,*) i,wv(i)
enddo
!read(5,*)
fmodbin=0.0d0 !initialize array
!resample with equal spacing
call binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
!write(0,*) "Done binning model"
deallocate(wmod,fmod) !get rid of uneven sampled grid
allocate(wmod(npt),fmod(npt)) !make new array with equal spaced grid
nmodel=npt
wmod=wv !copy work arrays
fmod=fmodbin
deallocate(wv,fmodbin) !get rid of work arrays

!lets fill in pixel values
xmax=xout*noversample !size of over sampled detector array.
ymax=yout*noversample
!array to hold detector array values
allocate(pixels(xmax,ymax),wpixels(xmax,ymax))
pixels=0.0d0 !initialize array to zero
dxmaxp1=dble(xmax+1) !xmax+1 converted to double
dymaxp1=dble(ymax+1) !ymax+1 converted to double
!allocate array for convolved image
allocate(cpixels(xmax,ymax),wcpixels(xmax,ymax))
cpixels=0.0d0 !initialize array to zero

do ntrace=1,ntracemax !loop over all traces
   write(0,*) "Trace # ",ntrace
   wpixels=0.0d0 !initialize array to zero
   !trace parts of array used for spectrum (speeds up convolution later)
   ybounds(1)=ymax !initalize to bad bounds for max/min search
   ybounds(2)=1
   !j=0
   do i=1,nmodel
      px=w2p(wmod(i),noversample,ntrace) !get x-pixel value
      !is the wavelength on the detector?
      if((px.gt.1.0d0).and.(px.lt.dxmaxp1))then
         py=ptrace(px,noversample,ntrace) !get y-pixel value
!      write(0,*) px,py
         if((py.gt.1.0d0).and.(py.lt.dymaxp1))then !check y-pixel value
            npx=int(px) !convert pixel values to integer
            npy=int(py)
            !find extremes of CCD use to speed up later
            ybounds(1)=min(ybounds(1),npy)
            ybounds(2)=max(ybounds(2),npy)
!           wpixels(npx,npy)=wpixels(npx,npy)+fmod(i) !add flux to pixel
!           add in response.

!           wmod is in A.. convert to nm for easy comparision
            awmod=wmod(i)/10.0d0
            if((awmod.lt.500.0).or.(awmod.gt.5500.0))then
               respond=0.0d0
            else
            ! cubic interpolation of response along traces
               select case(ntrace)
                  case(1)
                     call splint(ld,res1,yres1,nres,awmod,respond)
                  case(2)
                     call splint(ld,res2,yres2,nres,awmod,respond)
                  case(3)
                     call splint(ld,res3,yres3,nres,awmod,respond)
               end select
            endif
            fmodres=fmod(i)*respond !flux to add to pixel
            call addflux2pix(px,py,xmax,ymax,wpixels,fmodres)

!           j=j+1 !for counting number of wavelength samples
         endif
      endif
   enddo
   !write(0,*) "nsamples: ",j

   !Now we convolve the narrow spectra with the PSF kernel
   write(0,*) minval(wpixels),maxval(wpixels)
   wcpixels=0.0d0
   call convolve(xmax,ymax,wpixels,nrK,nKs,rKernel,noversample,wcpixels,&
      ybounds,ntrace)
   write(0,*) minval(wcpixels),maxval(wcpixels)

   !copy trace to master array for output
   pixels=pixels+wpixels !unconvolved image
   cpixels=cpixels+wcpixels !convolved image
enddo
deallocate(wcpixels,wpixels,yres1,yres2,yres3) !work arrays no longer needed

allocate(opixels(xout,yout)) !we resample the array for output
opixels=0.0d0 !initalize the array
dnossq=noversample*noversample
do i=noversample,xmax,noversample
   do j=noversample,ymax,noversample
      opixels(i/noversample,j/noversample)=                             &
         Sum(pixels(i-noversample+1:i,j-noversample+1:j))
!      opixels(i/noversample,j/noversample)=                             &
!         opixels(i/noversample,j/noversample)/dnossq
   enddo
enddo
fileout="spgen.fits"  !write out un-convolved 2D spectrum
call writefits(xout,yout,opixels,fileout) !make fits file.

opixels=0.0d0 !reinitialize the array
dnossq=noversample*noversample
do i=noversample,xmax,noversample  !resample (bin) the array.
   do j=noversample,ymax,noversample
      opixels(i/noversample,j/noversample)=                             &
         Sum(cpixels(i-noversample+1:i,j-noversample+1:j))
!      opixels(i/noversample,j/noversample)=                             &
!         opixels(i/noversample,j/noversample)/dnossq
   enddo
enddo
fileout="spgen_c.fits" !write out convolved 2D spectrum
call writefits(xout,yout,opixels,fileout) !make fits file.

end program spgen

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function ptrace(px,noversample,ntrace)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer noversample,i,ntrace
integer, parameter :: nc=5
real(double) :: ptrace,px,opx
real(double), dimension(nc) :: c1,c2,c3,c
data c1/275.685,0.0587943,-0.000109117,1.06605d-7,-3.87d-11/
data c2/254.109,-0.00121072,-1.84106e-05,4.81603e-09,-2.14646e-11/
data c3/203.104,-0.0483124,-4.79001e-05,0.0,0.0/

select case(ntrace)
   case(1)
      c=c1
   case(2)
      c=c2
   case(3)
      c=c3
end select

opx=px/dble(noversample) !account for oversampling
ptrace=c(1)
do i=2,nc
   !polynomial fit to trace. Good to about 1-2 pix
   ptrace=ptrace+opx**dble(i-1)*c(i)
enddo
ptrace=ptrace*dble(noversample) !account for oversampling

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function p2w(p,noversample,ntrace)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!convert pixel to wavelength
! n=1 polynomial fit
use precision
implicit none
integer :: i,noversample,ntrace
integer, parameter :: nc=5
real(double) :: pix,p,p2w
real(double), dimension(nc) :: c1,c2,c3,c
data c1/2.60188,-0.000984839,3.09333e-08,-4.19166e-11,1.66371e-14/
data c2/1.30816,-0.000480837,-5.21539e-09,8.11258e-12,5.77072e-16/
data c3/0.880545,-0.000311876,8.17443e-11,0.0,0.0/

select case(ntrace)
   case(1)
      c=c1
   case(2)
      c=c2
   case(3)
      c=c3
end select

pix=p/dble(noversample)
p2w=c(1)
do i=2,nc
   p2w=p2w+pix**dble(i-1)*c(i)
enddo
p2w=p2w*10000.0 !um->A

!write(0,*) p2w,pix,p
!read(5,*)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function w2p(w,noversample,ntrace)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!convert wavelength to pixel
! n=1 polynomial fit.
use precision
implicit none
integer :: i,noversample,ntrace
integer, parameter :: nc=5
real(double) w2p,w,wum
real(double), dimension(nc) :: c1,c2,c3,c
data c1/2957.38,-1678.19,526.903,-183.545,23.4633/
data c2/3040.35,-2891.28,682.155,-189.996,0.0/
data c3/2825.46,-3211.14,2.69446,0.0,0.0/

select case(ntrace)
   case(1)
      c=c1
   case(2)
      c=c2
   case(3)
      c=c3
end select

wum=w/10000.0 !A->um
w2p=c(1)
do i=2,nc
   w2p=w2p+wum**dble(i-1)*c(i)  !polynomial fit to trace. Good to about 1-2 pix
enddo
w2p=w2p*dble(noversample) !account for oversampling
!write(0,*) w2p,w/10000.0d0
!read(5,*)

return
end
