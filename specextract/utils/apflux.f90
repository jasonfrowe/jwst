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
real, allocatable, dimension(:) :: px,py

nap=20 !aperture for flux sums
nsky=50 !aperture for sky estimate

allocate(flux(naxes(1)))
k=1 !working on first trace for now.
do i=1,naxes(1)
   nm=max(1,int(dTrace(i,k))-nap)
   np=min(naxes(2),int(dTrace(i,k))+nap)
   flux(i)=Sum(Image(i,nm:np))
enddo

!plotting to have a look
allocate(px(naxes(1)),py(naxes(1)))
do i=1,naxes(1)
   px(i)=real(i)
enddo
py=real(flux)
!pmin=minval(py)
!py=log10(py-pmin+1.000)
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(minval(px),maxval(px),minval(py),maxval(py)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("X (pixels)","Y (Counts)","")
call pgline(naxes(1),px,py)
deallocate(px,py)

return
end subroutine apflux
