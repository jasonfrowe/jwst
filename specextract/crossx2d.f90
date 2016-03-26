program crossx2d
use precision
implicit none
integer iargc,nkeysmax,nxmax,nymax,nkeys,status,idcor,dumi
integer, dimension(2) :: naxes
real(double) :: Rmin,Rmax,bpix
real(double), allocatable, dimension(:,:) :: Image1,Image2,dark,tImage1,&
   tImage2
character(200) :: file1,file2,darkfile
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
interface
   subroutine xcorr2d(naxes,Image1,Image2)
      use precision
      implicit none
      integer, dimension(2) :: naxes
      real(double), dimension(:,:) :: Image1,Image2
   end subroutine xcorr2d
end interface

!Parameters
bpix=1.0e30  !mark bad pixels
nkeysmax=700 !maximum size of header
nxmax=2048   !maximum size of image
nymax=2048   !maximum size of image


if(iargc().lt.2)then
   write(0,*) "Usage: crossx2d <file1> <file2> [dark]"
   stop
endif

call getarg(1,file1)
call getarg(2,file2)
if(iargc().ge.3)then
   call getarg(3,darkfile)
   idcor=0
endif

!read in 1st FITS file
allocate(Image1(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image array failed.."
   write(0,*) "Status: ",status
   stop
endif
allocate(header(nkeysmax))
Image1=bpix !initialize Image array with bad-pixels
call getfits(file1,naxes,Image1,Rmin,Rmax,nkeys,header,bpix)
!write(0,*) "Rmin,Rmax: ",Rmin,Rmax

!read in 2nd FITS file
allocate(Image2(nxmax,nymax),stat=status)
if(status.gt.0) then !fix for gfortran
   write(0,*) "Allocation of Image array failed.."
   write(0,*) "Status: ",status
   stop
endif
Image2=bpix !initialize Image array with bad-pixels
call getfits(file2,naxes,Image2,Rmin,Rmax,nkeys,header,bpix)
!write(0,*) "Rmin,Rmax: ",Rmin,Rmax

if(idcor.eq.0)then !read in Dark and correct images if possible
   !read in dark FITS file
   allocate(dark(nxmax,nymax),stat=status)
   if(status.gt.0) then !fix for gfortran
      write(0,*) "Allocation of Image array failed.."
      write(0,*) "Status: ",status
      stop
   endif
   dark=bpix !initialize Image array with bad-pixels
   call getfits(darkfile,naxes,dark,Rmin,Rmax,nkeys,header,bpix)
!   write(0,*) "Rmin,Rmax: ",Rmin,Rmax

   !Dark correction
   Image1=Image1-dark
   Image2=Image2-dark

   !we don't need the dark anymore.
   deallocate(dark)
endif

!Rotate images to deal with CV3
if(naxes(1).lt.naxes(2))then
   Image1=transpose(Image1)
   Image2=transpose(Image2)
   dumi=naxes(1)
   naxes(1)=naxes(2)
   naxes(2)=dumi
endif

!resize images to minimize memory usage.
if((naxes(1).ne.nxmax).or.(naxes(2).ne.nymax))then
   allocate(tImage1(naxes(1),naxes(2)),tImage2(naxes(1),naxes(2)))
   tImage1(1:naxes(1),1:naxes(2))=Image1(1:naxes(1),1:naxes(2))
   tImage2(1:naxes(1),1:naxes(2))=Image2(1:naxes(1),1:naxes(2))
   deallocate(Image1,Image2)
   nxmax=naxes(1)
   nymax=naxes(2)
   allocate(Image1(nxmax,nymax),Image2(nxmax,nymax))
   Image1=tImage1
   Image2=tImage2
   deallocate(tImage1,tImage2)
endif

call xcorr2d(naxes,Image1,Image2)

end program crossx2d
