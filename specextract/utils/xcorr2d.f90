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
integer ( kind = 4 ) :: nh,nhs
real(double), allocatable, dimension(:,:) :: A,B,C,Ct
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: AC,BC,CC,ACs, &
   BCs
!local vars
integer :: noversample,XF,YF,i,j,imax,jmax,XFs,YFs,ii,jj
real(double) :: maxc


!choose oversampling and then allocate arrays
noversample=4  !64 needs 32G of memory!!
XF=naxes(1) !make arrays bigger for oversampling
YF=naxes(2)
allocate(A(XF,YF),B(XF,YF))
nh=(XF/2)+1
allocate(AC(nh,YF),BC(nh,YF))

A=Image1
B=Image2

!!**Insert spline here..

!FFT Plans for A&B
planA=fftw_plan_dft_r2c_2d(YF,XF,A,AC,FFTW_ESTIMATE)
planB=fftw_plan_dft_r2c_2d(YF,XF,B,BC,FFTW_ESTIMATE)

!execture FFT
call fftw_execute_dft_r2c(planA,A,AC)
call fftw_execute_dft_r2c(planB,B,BC)

!destroy plans
call fftw_destroy_plan(planA)
call fftw_destroy_plan(planB)

!deallocate A,B
deallocate(A,B)

!zero-pad AC,BC
XFs=naxes(1)*noversample
YFs=naxes(2)*noversample

!allocate oversampled ACs,BCs
nhs=XFs/2+1
allocate(ACs(nhs,YFs),BCs(nhs,YFs))
ACs=0.0d0 !initialize to zero for zero-padding
BCs=0.0d0

do i=1,(nh-1)/2
   do j=1,YF/2
      ACs(i,j)=AC(i,j)
      BCs(i,j)=BC(i,j)
   enddo
enddo
jj=YF*(noversample-1)
do i=1,(nh-1)/2
   do j=YF/2+1,YF
      ACs(i,j+jj)=AC(i,j)
      BCs(i,j+jj)=AC(i,j)
   enddo
enddo
ii=(nh-1)*(noversample-1)
do i=(nh-1)/2,nh
   do j=1,YF/2
      ACs(i+ii,j)=AC(i,j)
      BCs(i+ii,j)=BC(i,j)
   enddo
enddo
do i=(nh-1)/2,nh
   do j=YF/2+1,YF
      ACs(i+ii,j+jj)=AC(i,j)
      BCs(i+ii,j+jj)=BC(i,j)
   enddo
enddo

!open(unit=11,file="fft.dat")
!do i=1,nhs
!   write(11,'(4096(1PE17.10,1X))') (log10(abs(ACs(i,j))),j=1,YFs)
!enddo
!close(11)

!allocate C arrays
allocate(C(XFs,YFs),CC(nhs,YFs))

!calculate correlation
CC=ACs*dconjg(BCs)

!deallocate arrays
deallocate(ACs,BCs)

!FFT plan for C
planC=fftw_plan_dft_c2r_2d(YFs,XFs,CC,C,FFTW_ESTIMATE)
!initialize C
C=0.0d0
!invert FFT
call fftw_execute_dft_c2r(planC,CC,C)
!destroy plans
call fftw_destroy_plan(planC)
!deallocate arrays
deallocate(CC)
!simple normalization
C=C/maxval(C)

allocate(Ct(XFs,YFs))
ii=XFs/2
jj=YFs/2
do i=1,XFs/2
   do j=1,YFs/2
      Ct(i+ii,j+jj)=C(i,j)
   enddo
enddo
do i=XFs/2+1,XFs
   do j=1,YFs/2
      Ct(i-ii,j+jj)=C(i,j)
   enddo
enddo
do i=1,XFs/2
   do j=YFs/2+1,YFs
      Ct(i+ii,j-jj)=C(i,j)
   enddo
enddo
do i=XFs/2+1,XFs
   do j=YFs/2+1,YFs
      Ct(i-ii,j-jj)=C(i,j)
   enddo
enddo

!deallocate unneeded arrays
deallocate(C)

open(unit=11,file="fft.dat")
do i=1,XFs
   write(11,'(100000(1PE17.10,1X))') (Ct(i,j),j=1,YFs)
enddo
close(11)

maxc=0.0
do i=1,size(C(:,1))
   do j=1,size(C(1,:))
      if(Ct(i,j).gt.maxc)then
         maxc=Ct(i,j)
         imax=i
         jmax=j
      endif
   enddo
enddo
write(6,'(2(F9.5,1X))') dble(imax-XFs/2-1)/dble(noversample),           &
           dble(jmax-YFs/2-1)/dble(noversample)


end subroutine xcorr2d

!Old code
!real(double), allocatable, dimension(:) :: x1a,x2a
!real(double), allocatable, dimension(:,:) :: y2a
!real(double) :: x,y,f,dnover
!
!A=0.0d0 !initialize to zero (for zero-padding)
!B=0.0d0
!
!!variables for interpolation
!allocate(x1a(naxes(1)),x2a(naxes(2)),y2a(naxes(1),naxes(2)))
!do i=1,naxes(1) !create grid points based on indice
!   x1a(i)=dble(i)
!enddo
!do j=1,naxes(2) !create grid points based on indice
!   x2a(j)=dble(j)
!enddo
!!interpolate data for better cross-correlation
!write(0,*) "Start splie2"
!call splie2(x1a,x2a,Image1,naxes(1),naxes(2),y2a)
!
!write(0,*) "Start splin2"
!dnover=dble(noversample) !precompute int -> dble
!do i=1,XF
!   x=dble(i)/dnover  !interpolated grid
!   write(0,*) "i: ",i
!   do j=1,YF
!      y=dble(j)/dnover  !interpolated grid point
!      !find interpolated value
!      call splin2(x1a,x2a,Image1,y2a,naxes(1),naxes(2),x,y,f)
!      A(i,j)=f !construct interpolated matrix
!   enddo
!enddo
!
!!repeat splie2 and splin2 for Image2
!write(0,*) "Start splie2"
!call splie2(x1a,x2a,Image2,naxes(1),naxes(2),y2a)
!write(0,*) "Start splin2"
!do i=1,XF
!   x=dble(i)/dnover  !interpolated grid
!   do j=1,YF
!      y=dble(j)/dnover  !interpolated grid point
!      !find interpolated value
!      call splin2(x1a,x2a,Image2,y2a,naxes(1),naxes(2),x,y,f)
!      B(i,j)=f !construct interpolated matrix
!   enddo
!enddo
