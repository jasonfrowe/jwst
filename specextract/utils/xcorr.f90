subroutine xcorr(nlines,nTrace,apfluxl,apfluxu,avgsplitlength,tilt,imax)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
!iso_c_binding is for FFTW3 interface
use, intrinsic :: iso_c_binding
implicit none
!add in the FFTW3 modules
include 'fftw3.f03'

integer :: nlines,nTrace
real(double) :: avgsplitlength
real(double), dimension(:,:) :: apfluxl,apfluxu

integer :: XF,i,j,k,nplot,noversample,imax,jj
real, allocatable, dimension(:) :: px,py,bb
real(double) :: x,y,dnover,Pi,rad2deg,maxc,tilt,tilterr
real(double), allocatable, dimension(:) :: A,B,C,xl,yl,yu,yl2,yu2
character(80) :: tilttext

!FFTW3 vars
type(C_PTR) :: planA,planB,planC
integer ( kind = 4 ) :: nh
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: AC,BC,CC

!constants
Pi=acos(-1.d0)       !define Pi
rad2deg=180.0/Pi

k=1 !work on first trace

!choose oversampling and then allocate arrays
noversample=32
XF=nlines*noversample !make arrays bigger for oversampling
allocate(A(XF),B(XF),C(XF))
nh=(XF/2)+1
allocate(AC(nh),BC(nh),CC(nh))

A=0.0d0 !initialize to zero (for zero-padding)
B=0.0d0
!assign data to arrays for FFT
!A(1:nlines)=apfluxl(1:nlines,k)
!B(1:nlines)=apfluxu(1:nlines,k)

!interpolate data for better cross-correlation
allocate(xl(nlines),yl(nlines),yu(nlines),yl2(nlines),yu2(nlines))
do i=1,nlines
   xl(i)=dble(i)
   yl(i)=apfluxl(i,k)
enddo
call spline(xl,yl,nlines,1.d30,1.d30,yl2)
do i=1,nlines
   yu(i)=apfluxu(i,k)
enddo
call spline(xl,yu,nlines,1.d30,1.d30,yu2)

dnover=dble(noversample)
do i=1,XF
!   write(0,*) "i: ",i
   x=dble(i)/dnover
   call splint(xl,yl,yl2,nlines,x+0.0,y)
!   write(0,*) "y: ",y
   A(i)=y
   call splint(xl,yu,yu2,nlines,x+0.0,y)
   B(i)=y
enddo

!generate FFTW plans
planA=fftw_plan_dft_r2c_1d(XF,A,AC,FFTW_ESTIMATE)
planB=fftw_plan_dft_r2c_1d(XF,B,BC,FFTW_ESTIMATE)
planC=fftw_plan_dft_c2r_1d(XF,CC,C,FFTW_ESTIMATE)

!execute FFT
call fftw_execute_dft_r2c(planA,A,AC)
call fftw_execute_dft_r2c(planB,B,BC)
!destroy plans
call fftw_destroy_plan(planA)
call fftw_destroy_plan(planB)
!calculate correlation
CC=AC*dconjg(BC)
!invert FFT
call fftw_execute_dft_c2r(planC,CC,C)
!C=C/dble(nlines)
C=C/maxval(C)

maxc=0.0
do i=1,size(C)
   if(C(i).gt.maxc)then
      maxc=C(i)
      imax=i
   endif
enddo


jj=size(C)
if(imax.gt.jj/2)then
   tilt=atan(real(imax-jj)/dnover/avgsplitlength)*rad2deg
   write(0,*) "pixels: ",dble(imax-jj)/dnover
else
   tilt=atan(real(imax)/dnover/avgsplitlength)*rad2deg
   write(0,*) "pixels: ",dble(imax)/dnover
endif
tilterr=atan(-1.0d0/dnover/avgsplitlength)*rad2deg
write(0,*) "tilt,err: ",tilt,tilterr
!tilt=-tilt
!write(6,'(F9.4,1X,I5,1X,A12)') tilt,imax,"(deg,pixels)"
!write(0,*) tilt,imax,size(C)
if(imax.gt.size(C)/2)then
   imax=imax-size(C)
endif
write(tilttext,'(A6,F9.4)') "Tilt: ",tilt

!j=50*noversample
j=min(50,nlines)*noversample
allocate(px(j*2),py(j*2),bb(4))
nplot=0
k=1 !plot first Trace
do i=size(C)-j+1,size(C)
   nplot=nplot+1
!   px(nplot)=real(nplot-j)
   px(nplot)=atan(real(nplot-j)/dnover/avgsplitlength)*rad2deg
   py(nplot)=real(C(i))
enddo
do i=1,j
   nplot=nplot+1
!   px(nplot)=real(i)
   px(nplot)=atan(real(i)/dnover/avgsplitlength)*rad2deg
   py(nplot)=real(C(i))
enddo
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.20,0.95) !make room around the edges for labels
bb(1)=minval(px(1:nplot))
bb(2)=maxval(px(1:nplot))
bb(3)=minval(py(1:nplot))
bb(4)=maxval(py(1:nplot))
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
!call pgwindow(50.0,250.0,minval(py),maxval(py(50:250)))
call pgsch(2.0)
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("Tilt (deg)","Corr","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "Tilt (deg)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Corr")
call pgline(nplot,px,py)
call pgwindow(0.0,1.0,0.0,1.0)
call PGPTXT (0.8, 0.8, 0.0, 0.0, tilttext)
deallocate(px,py)

return
end subroutine xcorr
