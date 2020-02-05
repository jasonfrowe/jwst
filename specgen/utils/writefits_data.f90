subroutine writefitsdata(funit,xout,yout,pixels,ngroup,nint)
use precision
implicit none
!import arrays
integer :: funit,xout,yout,ngroup,nint
real(double), dimension(:,:) :: pixels
!local arrays
integer :: i,j,k,l !counters
integer :: status,bitpix,naxis,npixels,group,firstpix,nbuf,nbuffer
integer, dimension(:), allocatable :: naxes,buffer
real(double) :: dngrpfac
!pixels for writing
real(double), dimension(:,:,:,:), allocatable :: pixelsout
!plot arrrays
real(double) :: bpix,tavg,sigscale

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
!display fits file
!call pgopen('?')
call pgopen('/xserve')
!call pgopen('trace.ps/vcps')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
bpix=1.0e30
tavg=0.0
sigscale=3.0
call pgpage()
call displayfits(xout,yout,pixels,bpix,tavg,sigscale)
call pgclos()


status=0 !tracks errors for FITSIO routines

!BITPIX = 16 means that the image pixels will consist of 16-bit
!integers.  The size of the image is given by the NAXES values.
bitpix=-32 !using float-doubles. 
naxis=4 !JWST obs have 4 axes
allocate(naxes(naxis)) 
naxes(1) = xout
naxes(2) = yout
naxes(3) = ngroup
naxes(4) = nint


allocate(pixelsout(xout,yout,ngroup,nint))
do k=1,ngroup
   dngrpfac=dble(ngroup-k+1) 
   pixelsout(1:xout,1:yout,k,1)=pixels(1:xout,1:yout)/dngrpfac
enddo

!insert a new IMAGE extension immediately following the CHDU
call FTIIMG(funit,bitpix,naxis,naxes,status)

firstpix=1
group=1 !this var does nothing, leave it alone
npixels=1
do i=1,naxis
	npixels=npixels*naxes(i)
enddo
call ftpprd(funit,group,firstpix,npixels,pixelsout,status)

return
end subroutine writefitsdata
