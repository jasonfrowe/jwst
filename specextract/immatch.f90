program immatch
use precision
implicit none
integer :: iargc,nkeysmax,nxmax,nymax,status,dumi,nkeys,i,j
integer, dimension(2) :: naxes1,naxes2
real(double) :: bpix,Rmin,Rmax,tavg
real(double), allocatable, dimension(:) :: a
real(double), allocatable, dimension(:,:) :: Image1,Image2,tImage,Image3
character(80) :: file1,file2,fileout
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
   subroutine writefits2(nxmax,nymax,parray,bpix,tavg,nkeys,header,fileout)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax,nkeys
      real(double), intent(inout) :: bpix,tavg
      real(double), dimension(:,:), intent(inout) :: parray
      character(80), intent(inout) :: fileout
      character(80), dimension(:), intent(inout) :: header
   end subroutine writefits2
end interface

!get filename
if(iargc().lt.2)then
   write(0,*) "Usage: immatch <FITS1> <FITS2>"
   write(0,*) "  <FITS1>   : FITS file containing SOSS data"
   write(0,*) "  <FITS2>   : FITS file containing SOSS data"
   write(0,*) " Program calculates FITS1-scale*(FITS2+ZPT)"
   stop
endif
call getarg(1,file1)
call getarg(2,file2)

!read in FITS file
bpix=1.0e30  !mark bad pixels
nkeysmax=700
nxmax=2048
nymax=2048
allocate(Image1(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image1 array failed.."
   write(0,*) "Status: ",status
   stop
endif
allocate(header(nkeysmax))
Image1=bpix !initialize Image array with bad-pixels
call getfits(file1,naxes1,Image1,Rmin,Rmax,nkeys,header,bpix)
!write(0,*) "FITS read..",naxes
!For CV3 we need to transpose the Image.
if(naxes1(1).lt.naxes1(2))then
   Image1=transpose(Image1)
   dumi=naxes1(1)
   naxes1(1)=naxes1(2)
   naxes1(2)=dumi
endif

!compact memory requirements for array
if((naxes1(1).ne.nxmax).or.(naxes1(2).ne.nymax))then
   allocate(tImage(naxes1(1),naxes1(2)))
   tImage(1:naxes1(1),1:naxes1(2))=Image1(1:naxes1(1),1:naxes1(2))
   deallocate(Image1)
   allocate(Image1(naxes1(1),naxes1(2)))
   Image1=tImage
   deallocate(tImage)
endif

!read in second Image
allocate(Image2(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image2 array failed.."
   write(0,*) "Status: ",status
   stop
endif
Image2=bpix !initialize Image array with bad-pixels
call getfits(file2,naxes2,Image2,Rmin,Rmax,nkeys,header,bpix)
!write(0,*) "FITS read..",naxes
!For CV3 we need to transpose the Image.
if(naxes2(1).lt.naxes2(2))then
   Image2=transpose(Image2)
   dumi=naxes2(1)
   naxes2(1)=naxes2(2)
   naxes2(2)=dumi
endif

!compact memory requirements for array
if((naxes2(1).ne.nxmax).or.(naxes2(2).ne.nymax))then
   allocate(tImage(naxes2(1),naxes2(2)))
   tImage(1:naxes2(1),1:naxes2(2))=Image2(1:naxes2(1),1:naxes2(2))
   deallocate(Image2)
   allocate(Image2(naxes2(1),naxes2(2)))
   Image2=tImage
   deallocate(tImage)
endif

!check for dimension match
if((naxes1(1).ne.naxes2(1)).or.(naxes1(2).ne.naxes2(2)))then
   write(0,*) "Error: Dimension of Images does not match"
   write(0,*) "Image1: ",naxes1(1),naxes1(2)
   write(0,*) "Image2: ",naxes2(1),naxes2(2)
endif

!update physical bounds of Image arrays
nxmax=naxes1(1)
nymax=naxes1(2)

!display fits file
call pgopen('?')
!call pgopen('/xserve')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage()
tavg=0.0 !displays a time on the image
call displayfits(nxmax,nymax,Image1,bpix,tavg,0.0d0)
call PGSFS(2)
CALL PGSCR(16, 1.0, 1.0, 1.0)
call pgsci(16)
call pgrect(1675.0,2038.0,140.0,220.0)
call pgsci(1)
call pgpage()
call displayfits(nxmax,nymax,Image2,bpix,tavg,0.0d0)
CALL PGSCR(16, 1.0, 1.0, 1.0)
call pgsci(16)
call pgrect(1675.0,2038.0,140.0,220.0)
call pgsci(1)
call PGSFS(1)

!open(unit=11,file="match.dat")
!do i=1675,2038
!   do j=140,220
!      write(11,*) Image1(i,j),Image2(i,j)
!   enddo
!enddo
!close(11)

allocate(a(2))
a(1)=50.4926
a(2)=0.722905

allocate(Image3(nxmax,nymax))
Image3=Image1-(Image2*a(2)+a(1))

call pgpage()
call displayfits(nxmax,nymax,Image3,bpix,tavg,0.0d0)

fileout="test.fits"
call writefits2(nxmax,nymax,Image3,bpix,tavg,nkeys,header,fileout)

call pgclos()

end program immatch
