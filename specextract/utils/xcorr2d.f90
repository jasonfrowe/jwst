subroutine xcorr2d(naxes,Image1,Image2)
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'
!import vars
integer, dimension(2) :: naxes
real(double), dimension(:,:) :: Image1,Image2
!FFTW3 vars
type(C_PTR) :: planA,planB,planC
integer ( kind = 4 ) :: nh
real(double), allocatable, dimension(:,:) :: A,B,C
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: AC,BC,CC
!local vars
integer :: noversample,XF,YF,i,j,imax,jmax
real(double) :: dnover,x,y,f,maxc
real(double), allocatable, dimension(:) :: x1a,x2a
real(double), allocatable, dimension(:,:) :: y2a

!choose oversampling and then allocate arrays
noversample=2
XF=naxes(1)*noversample !make arrays bigger for oversampling
YF=naxes(2)*noversample
allocate(A(XF,YF),B(XF,YF),C(XF,YF))
nh=(XF/2)+1
allocate(AC(nh,YF),BC(nh,YF),CC(nh,YF))

A=0.0d0 !initialize to zero (for zero-padding)
B=0.0d0

!variables for interpolation
allocate(x1a(naxes(1)),x2a(naxes(2)),y2a(naxes(1),naxes(2)))
do i=1,naxes(1) !create grid points based on indice
   x1a(i)=dble(i)
enddo
do j=1,naxes(2) !create grid points based on indice
   x2a(j)=dble(j)
enddo

!interpolate data for better cross-correlation
write(0,*) "Start splie2"
call splie2(x1a,x2a,Image1,naxes(1),naxes(2),y2a)

write(0,*) "Start splin2"
dnover=dble(noversample) !precompute int -> dble
do i=1,XF
   x=dble(i)/dnover  !interpolated grid
   write(0,*) "i: ",i
   do j=1,YF
      y=dble(j)/dnover  !interpolated grid point
      !find interpolated value
      call splin2(x1a,x2a,Image1,y2a,naxes(1),naxes(2),x,y,f)
      A(i,j)=f !construct interpolated matrix
   enddo
enddo

!repeat splie2 and splin2 for Image2
write(0,*) "Start splie2"
call splie2(x1a,x2a,Image2,naxes(1),naxes(2),y2a)
write(0,*) "Start splin2"
do i=1,XF
   x=dble(i)/dnover  !interpolated grid
   do j=1,YF
      y=dble(j)/dnover  !interpolated grid point
      !find interpolated value
      call splin2(x1a,x2a,Image2,y2a,naxes(1),naxes(2),x,y,f)
      B(i,j)=f !construct interpolated matrix
   enddo
enddo

!FFT Plans
planA=fftw_plan_dft_r2c_2d(YF,XF,A,AC,FFTW_ESTIMATE)
planB=fftw_plan_dft_r2c_2d(YF,XF,B,BC,FFTW_ESTIMATE)
planC=fftw_plan_dft_c2r_2d(YF,XF,CC,C,FFTW_ESTIMATE)

!execture FFT
call fftw_execute_dft_r2c(planA,A,AC)
call fftw_execute_dft_r2c(planB,B,BC)

!destroy plans
call fftw_destroy_plan(planA)
call fftw_destroy_plan(planB)

!calculate correlation
CC=AC*dconjg(BC)

!invert FFT
call fftw_execute_dft_c2r(planC,CC,C)
!destroy plans
call fftw_destroy_plan(planC)

maxc=0.0
do i=1,size(C(:,1))
   do j=1,size(C(1,:))
      if(C(i,j).gt.maxc)then
         maxc=C(i,j)
         imax=i
         jmax=j
      endif
   enddo
enddo
write(0,*) "shift: ",maxc,imax,jmax

end subroutine xcorr2d
