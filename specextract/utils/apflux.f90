subroutine apflux(naxes,Image,bpix,nTrace,dTrace)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nTrace,i,nap,nsky,k,nm,np
integer, dimension(2) :: naxes
real(double) :: bpix
real(double), dimension(:,:) :: Image,dTrace
real(double), allocatable, dimension(:) :: flux
!plotting
real :: pmin
real, allocatable, dimension(:) :: px,py,bb

nap=20 !aperture for flux sums
nsky=50 !aperture for sky estimate

allocate(flux(naxes(1)))
k=1 !working on first trace for now.
do i=1,naxes(1)
   nm=max(1,int(dTrace(i,k))-nap)
   np=min(naxes(2),int(dTrace(i,k))+nap)
   flux(i)=Sum(Image(i,nm:np))
enddo

open(unit=11,file='apflux.dat')
do i=1,naxes(1)
   write(11,*) i,flux(i)
enddo
close(11)

!plotting to have a look
allocate(px(naxes(1)),py(naxes(1)),bb(4))
do i=1,naxes(1)
   px(i)=real(i)
enddo
py=real(flux)
!pmin=minval(py)
!py=log10(py-pmin+1.000)
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
bb(1)=minval(px)
bb(2)=maxval(px)
bb(3)=minval(py)
bb(4)=maxval(py)
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "X (pixels)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Y (counts)")
call pgline(naxes(1),px,py)
deallocate(px,py)

return
end subroutine apflux
