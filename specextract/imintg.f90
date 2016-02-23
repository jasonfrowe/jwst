program imintg
!reads in data cubes and integrates.
use precision
implicit none
integer :: iargc,i
integer, allocatable, dimension(:) :: naxes,imax
real(double) :: bpix
real(double), allocatable, dimension(:,:,:) :: Imagecube
character(80) :: filename

interface
   subroutine readcfits(filename,naxes,imax,Imagecube,bpix)
      use precision
      implicit none
      integer, dimension(:) :: naxes,imax
      real(double) :: bpix
      real(double), dimension(:,:,:) :: Imagecube
      character(80) :: filename
   end subroutine
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: imintg"
   stop
endif

call getarg(1,filename)

allocate(naxes(3),imax(3))
imax(1)=2048
imax(2)=2048
imax(3)=100
allocate(Imagecube(imax(1),imax(2),imax(3)))
bpix=-1.0e10 !marking bad pixels
call readcfits(filename,naxes,imax,Imagecube,bpix)

end program imintg
