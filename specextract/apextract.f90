program apextract
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: iargc,nkeysmax,nxmax,nymax,status,nkeys
integer, dimension(2) :: naxes
real(double) :: bpix,Rmin,Rmax
real(double), allocatable, dimension(:,:) :: Image
character(80) :: Imagename,tracename
character(80), allocatable, dimension(:) :: header

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

!get filename
if(iargc().lt.2)then
   write(0,*) "Usage: apextract <FITS> <tracepsf.dat>"
   write(0,*) "  <FITS>         : FITS file containing SOSS data"
   write(0,*) "  <tracepsf.dat> : Output from spectextract"
   stop
endif
call getarg(1,Imagename)
call getarg(2,tracename)

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


end program apextract
