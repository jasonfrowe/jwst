subroutine convolveft(xmax,ymax,pixels,nrK,nKs,rKernel,noversample,       &
   cpixels,ybounds,ntrace)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'

integer :: xmax,ymax,nrK,nKs,noversample,ntrace,XF,YF,i
integer, dimension(2) :: ybounds
real(double) :: wl,wls,wle,dwl,p,wlp,p2w,fac
real(double), dimension(xmax,ymax) :: pixels,cpixels
real(double), dimension(nrK,nKs,nKs) :: rKernel
real(double), allocatable, dimension(:,:) :: Kernel,A,B,C
character(80) :: line

!FFTW3 vars
type(C_PTR) :: plan
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: AC,BC,CC

interface
   subroutine genkernel(nrK,nKs,rKernel,Kernel,wl,wls,wle,dwl)
      use precision
      integer, intent(inout) :: nrK,nKs
      real(double), intent(inout) :: wl,wls,wle,dwl
      real(double), dimension(:,:),intent(inout) :: Kernel
      real(double), dimension(:,:,:), intent(inout) :: rKernel
   end subroutine
end interface

cpixels=0.0d0! initalize convolved image to zero.

wls=0.5d0 !starting wavelength of Kernels
wle=3.4d0 !ending wavelengths of Kernels
dwl=0.1d0 !wavelength intervals
allocate(Kernel(nKs,nKs))
XF=nKs+xmax !size of zero padded arrays for convolution
YF=nKs+ymax
allocate(A(XF,YF),B(XF,YF),C(XF,YF)) !arrays to apply FFT
nh=(XF/2)+1 !size for complex array
allocate(AC(nh,YF),BC(nh,YF),CC(nh,YF)) !allocate complex arrays for FT

wl=wls
do while(wl.lt.wle)
   write(line,502) "Convolution: ",wl,wle
   502 format(A13,F6.3,1X,F6.3)
   call ovrwrt(line,2)
!   write(0,*) "wl: ",wl
!gets get a Kernel (wl is wavelength in um)
   call genkernel(nrK,nKs,rKernel,Kernel,wl,wls,wle,dwl)
   Kernel=transpose(Kernel)

   A=0.0d0 !initalize to zero
   B=0.0d0
   C=0.0d0
   A(1:xmax,1:ymax)=pixels(1:xmax,1:ymax) !assign Image to zero-padded A
!B(1:nKs/2,1:nKs/2)=Kernel(nKs/2:nKs,nKs/2:nKs) !assign Kernel to zero-padded B
!B(nKs/2:nKs,nKs/2:nKs)=Kernel(1:nKs/2,1:nKs/2)
   B(1:nKs,1:nKs)=Kernel(1:nKs,1:nKs)
!   write(0,*) "Computing FFTs"
   plan=fftw_plan_dft_r2c_2d(YF,XF,A,AC,FFTW_ESTIMATE)
   call fftw_execute_dft_r2c(plan,A,AC)
   call fftw_destroy_plan(plan)
   plan=fftw_plan_dft_r2c_2d(YF,XF,B,BC,FFTW_ESTIMATE)
   call fftw_execute_dft_r2c(plan,B,BC)
   call fftw_destroy_plan(plan)
   !multiply
   CC=AC*BC
!   write(0,*) "Start iFFT"
   plan=fftw_plan_dft_c2r_2d(YF,XF,CC,C,FFTW_ESTIMATE)
   call fftw_execute_dft_c2r(plan,CC,C)
   call fftw_destroy_plan(plan)
   do i=1,xmax
      p=dble(i)
      wlp=p2w(p,noversample,ntrace)/10000.0d0
      if(abs(wlp-wl).gt.dwl)then
         cpixels(i,1:ymax)=cpixels(i,1:ymax)+0.0d0
      else
         fac=1.0d0-abs(wlp-wl)/dwl
         cpixels(i,1:ymax)=cpixels(i,1:ymax)+C(i,1:ymax)/dble(xmax*ymax)*fac
      endif
   enddo
   wl=wl+dwl
enddo


return
end subroutine convolveft










