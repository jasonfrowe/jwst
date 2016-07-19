program dispfits
use precision
implicit none
integer :: iargc,i,nkeys,nkeysmax,nxmax,nymax,status,dumi
integer, dimension(2) :: naxes
real(double) :: bpix,Rmin,Rmax,tavg,sigscale
real(double), allocatable, dimension(:) :: sigs
real(double), allocatable, dimension(:,:) :: Image,tImage
character(80) :: filename
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

if(iargc().lt.1)then
   write(0,*) "Usage: dispfits <Image1> [Image2] ... [ImageN]"
   write(0,*) "  <Image1> : First FITS image to display"
   write(0,*) "  [Image2] : Second FITS image to display (optional)"
   write(0,*) "  [ImageN] : Nth FITS image to display (optional)"
   stop
endif

!display fits file
call pgopen('?')
!call pgopen('/xserve')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)

nkeysmax=700
allocate(header(nkeysmax),sigs(iargc()))
sigs=1000.0d0
sigs(2)=20.0d0

do i=1,iargc()  !loop over command-line arguments

   call getarg(i,filename) !get filename
   write(0,'(A80)') filename
   !read in FITS file
   bpix=1.0e30  !mark bad pixels

   nxmax=2048
   nymax=2048
   allocate(Image(nxmax,nymax),stat=status)
   if(status.gt.0) then
      write(0,*) "Allocation of Image1 array failed.."
      write(0,*) "Status: ",status
      stop
   endif

   Image=bpix !initialize Image array with bad-pixels
   call getfits(filename,naxes,Image,Rmin,Rmax,nkeys,header,bpix)
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
      allocate(Image(naxes(1),naxes(2)))
      Image=tImage
      deallocate(tImage)
   endif
   nxmax=naxes(1)
   nymax=naxes(2)

   call pgpage()
   tavg=0.0 !displays a time on the image
   sigscale=sigs(i)
   call displayfits(nxmax,nymax,Image,bpix,tavg,sigscale)

   deallocate(Image)

enddo


call pgclos()

end program dispfits
