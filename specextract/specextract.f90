program specextract
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nkeys,nkeysmax,nxmax,nymax,iargc,i,nline,nTrace,j,status,    &
   nfit,dumi,nap,nsky
integer, dimension(2) :: naxes
real :: pmin
integer :: nplot
real, allocatable, dimension(:) :: px,py,bb
real(double) :: Rmin,Rmax,bpix,tavg
real(double), allocatable, dimension(:) :: posguess,flux
real(double), allocatable, dimension(:,:) :: Image,tImage,dTrace,bf,solpsf
character(80) :: Imagename,cline
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
   subroutine displayfits(nxmax,nymax,parray,bpix,tavg,sigscale)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray
      real(double), intent(inout) :: bpix,tavg
      real(double), intent(in) :: sigscale
   end subroutine displayfits
end interface
interface
   subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf,solpsf,posguess)
      use precision
      implicit none
      integer :: nline,nTrace
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: bpix
      real(double), dimension(:) :: posguess
      real(double), dimension(:,:), intent(inout) :: Image,dTrace,bf,solpsf
   end subroutine trace
end interface
interface
   subroutine apflux(naxes,Image,bpix,nTrace,dTrace,nap,nsky,flux)
      use precision
      implicit none
      integer, intent(inout) :: nTrace,nap,nsky
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: bpix
      real(double), dimension(:), intent(inout) :: flux
      real(double), dimension(:,:), intent(inout) :: Image,dTrace
   end subroutine apflux
end interface

!parameters
nline=800  !magic line for trace
nTrace=3   !number of traces to track

!get filename
if(iargc().lt.2)then
   write(0,*) "Usage: specextract <FITS> <ntrace>"
   write(0,*) "  <FITS>   : FITS file containing SOSS data"
   write(0,*) "  <ntrace> : number of orders to trace.  Must be 1 or greater"
   write(0,*) "  [nStart] : column to start trace with, optional (default=800)"
   stop
endif
call getarg(1,Imagename)
call getarg(2,cline)
read(cline,*) ntrace  !read in ntrace
if(ntrace.lt.1)then   !check that value is valid
   write(0,*) "<ntrace> must be greater than zero"
   stop
endif
write(0,*) "nTrace: ",ntrace

!read in FITS file
bpix=1.0e30  !mark bad pixels
nkeysmax=700
nxmax=2048
nymax=2048
allocate(Image(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image array failed.."
   write(0,*) "Status: ",status
   stop
endif
allocate(header(nkeysmax))
Image=bpix !initialize Image array with bad-pixels
call getfits(Imagename,naxes,Image,Rmin,Rmax,nkeys,header,bpix)
!write(0,*) "FITS read..",naxes
!For CV3 we need to transpose the Image.
if(naxes(1).lt.naxes(2))then
   Image=transpose(Image)
   dumi=naxes(1)
   naxes(1)=naxes(2)
   naxes(2)=dumi
endif

!compact memory requirements for array
if((naxes(1).ne.nxmax).or.(naxes(2).ne.nymax))then
   allocate(tImage(naxes(1),naxes(2)))
   tImage(1:naxes(1),1:naxes(2))=Image(1:naxes(1),1:naxes(2))
   deallocate(Image)
   nxmax=naxes(1)
   nymax=naxes(2)
   allocate(Image(nxmax,nymax))
   Image=tImage
   deallocate(tImage)
endif

!read in magic line if provided and check that value is valid
if(iargc().ge.3)then
   call getarg(3,cline)
   read(cline,*) nline
endif
if((nline.lt.1).or.(nline.gt.naxes(1)))then
   write(0,'(A30,I5)') "nline invalid, 1 <= nline <= ",naxes(1)
   stop
endif

!if there is a magic line, then one can also allow initial guesses
!for trace location.  Good for low S/N cases
allocate(posguess(ntrace))
posguess=-1.0d0 !initiate to negative number to indicate invalid value
do i=1,ntrace
   if(iargc().ge.3+i)then
      call getarg(3+i,cline)
      read(cline,*) posguess(i)
      if((posguess(i).le.0.0).or.(posguess(i).gt.naxes(2)))then
         write(0,'(A9,I2,A31)') "position ",i,"invalid, will be safely ignored"
      endif
   endif
enddo

!no need remove negative pixels
!do i=1,naxes(1)
!   do j=1,naxes(2)
!      Image(i,j)=max(-0.0d0,Image(i,j))
!   enddo
!enddo

!write(0,*) "Naxes:" ,naxes(1),naxes(2)
!write(0,*) "Rmin, Rmax: ",Rmin,Rmax

!display fits file
!call pgopen('?')
call pgopen('/xserve')
!call pgopen('trace.ps/vcps')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage()
tavg=0.0 !displays a time on the image
call displayfits(nxmax,nymax,Image,bpix,tavg,20.0d0)

!plot sum of each column
allocate(px(nxmax),py(nxmax),bb(4))
do i=1,nxmax
   px(i)=real(i)
enddo
py=real(Sum(Image,2))
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
bb(1)=0.0
bb(2)=real(nxmax)
bb(3)=minval(py)
bb(4)=maxval(py)
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "X (pixels)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Y (counts)")
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
bb(1)=170.0 !minval(px)
bb(2)=250.0 !maxval(px)
bb(3)=minval(py)
bb(4)=maxval(py(120:250))
call pgwindow(bb(1),bb(2),bb(3),bb(4))
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "X (pixels)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Y (counts)")
call pgline(nymax,px,py)
deallocate(px,py)

!generate trace
nfit=1+9*ntrace
allocate(dtrace(naxes(1),nTrace),bf(naxes(1),nTrace),solpsf(naxes(1),nfit))
dTrace=0.0d0
call trace(naxes,Image,bpix,nline,nTrace,dTrace,bf,solpsf,posguess)
write(6,'(I2,1X,I4)') ntrace,naxes(1)
do i=1,naxes(1)
!   write(6,'(I4,3(1X,F11.3),3(1X,1PE17.10))') i,(dTrace(i,j),j=1,3),    &
!      (bf(i,j),j=1,3)
    write(6,'(I4,90(1X,1PE17.10))') i,(solpsf(i,j),j=1,nfit)
enddo

!extract aperture
nap=20 !aperture for flux sums
nsky=50 !aperture for sky estimate
allocate(flux(naxes(1)))
call apflux(naxes,Image,bpix,nTrace,dTrace,nap,nsky,flux)

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

!plot Image again
call pgpanl(1,1)
call displayfits(nxmax,nymax,Image,bpix,tavg,20.0d0)
allocate(px(naxes(1)),py(naxes(1)))
do i=1,nTrace
   nplot=0
   do j=1,naxes(1)
      if(bf(j,i).gt.1.0)then
         nplot=nplot+1
         px(nplot)=real(j)
         py(nplot)=real(dTrace(j,i))
      endif
   enddo
   call pgsci(2+i)
   call pgline(nplot,px,py)
enddo


call pgclos()

end program specextract
