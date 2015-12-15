program specextract
use precision
implicit none
integer :: nkeys,nkeysmax,nxmax,nymax,iargc,i,nline,nTrace,j,status
integer, dimension(2) :: naxes
real :: pmin
integer :: nplot
real, allocatable, dimension(:) :: px,py
real(double) :: Rmin,Rmax,bpix,tavg
real(double), allocatable, dimension(:,:) :: Image,tImage,dTrace,bf
character(80) :: Imagename
character(80), allocatable, dimension(:) :: header

!Interfaces to subroutines
interface
   subroutine getfits(Refname,naxes,Ref,Rmin,Rmax,nkeys,header,bpix)
      use precision
      implicit none
      integer :: nkeys
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: Rmin,Rmax,bpix
      real(double), dimension(:,:), intent(inout) :: Ref
      character(80), intent(inout) :: Refname
      character(80), dimension(:), intent(inout) :: header
   end subroutine getfits
end interface
interface
   subroutine displayfits(nxmax,nymax,parray,bpix,tavg)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray
      real(double), intent(inout) :: bpix,tavg
   end subroutine displayfits
end interface
interface
   subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf)
      use precision
      implicit none
      integer :: nline,nTrace
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: bpix
      real(double), dimension(:,:), intent(inout) :: Image,dTrace,bf
   end subroutine trace
end interface
interface
   subroutine apflux(naxes,Image,bpix,nTrace,dTrace)
      use precision
      implicit none
      integer, intent(inout) :: nTrace
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: bpix
      real(double), dimension(:,:), intent(inout) :: Image,dTrace
   end subroutine apflux
end interface

!parameters
nline=800  !magic line for trace
nTrace=3   !number of traces to track


!get filename
if(iargc().lt.1)then
   write(0,*) "Usage: specextract FITS"
   write(0,*) "  FITS : FITS file containing SOSS data"
   stop
endif
call getarg(1,Imagename)

!read in FITS file
bpix=1.0e30  !mark bad pixels
nkeysmax=700
nxmax=2048
nymax=512
allocate(Image(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image array failed.."
   write(0,*) "Status: ",status
   stop
endif
allocate(header(nkeysmax))
Image=bpix !initialize Image array with bad-pixels
call getfits(Imagename,naxes,Image,Rmin,Rmax,nkeys,header,bpix)

do i=1,naxes(1)
   do j=1,naxes(2)
      Image(i,j)=max(-0.0d0,Image(i,j))
   enddo
enddo

!write(0,*) "Naxes:" ,naxes(1),naxes(2)
!write(0,*) "Rmin, Rmax: ",Rmin,Rmax

!display fits file
call pgopen('?')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage()
tavg=0.0 !displays a time on the image
call displayfits(nxmax,nymax,Image,bpix,tavg)

!plot each column
allocate(px(nxmax),py(nxmax))
do i=1,nxmax
   px(i)=real(i)
enddo
py=real(Sum(Image,2))
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(0.0,real(nxmax),minval(py),maxval(py)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("X (pixels)","Y (Counts)","")
call pgline(nxmax,px,py)
deallocate(px,py)

allocate(px(nymax),py(nymax))
do i=1,nymax
   px(i)=real(i)
enddo
py=real(Image(nline,:))
!pmin=minval(py)
!py=log10(py-pmin+1.0d0)
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
!call pgwindow(minval(px),maxval(px),minval(py),maxval(py)) !plot scale
call pgwindow(200.0,350.0,minval(py),maxval(py))
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("X (pixels)","Y (Counts)","")
call pgline(nymax,px,py)
deallocate(px,py)

!generate trace
allocate(dtrace(naxes(1),nTrace),bf(naxes(1),nTrace))
dTrace=0.0d0
call trace(naxes,Image,bpix,nline,nTrace,dTrace,bf)

!extract aperture
!nTrace=1
call apflux(naxes,Image,bpix,nTrace,dTrace)
!nTrace=3

!plot Image again
call pgpanl(1,1)
call displayfits(nxmax,nymax,Image,bpix,tavg)
allocate(px(naxes(1)),py(naxes(1)))
do i=1,3!nTrace
   nplot=0
   do j=1,naxes(1)
      if(bf(j,i).gt.1000.0)then
         nplot=nplot+1
         px(nplot)=real(j)
         py(nplot)=real(dTrace(nplot,i))
      endif
   enddo
   call pgsci(2+i)
   call pgline(nplot,px,py)
enddo


call pgclos()

end program specextract
