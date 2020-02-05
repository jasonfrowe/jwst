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

!insert a new IMAGE extension immediately following the CHDU
call FTIIMG(funit,bitpix,naxis,naxes,status)

!Write the array to the FITS file.
allocate(buffer(yout))
firstpix=1
group=1 !this var does nothing, leave it alone
l=1 !this controls nint (if we have multiple images or resets)
do k=1,ngroup
   dngrpfac=dble(ngroup-k+1) !we are implementing a strick linear ramp for this simulation

   npixels=naxes(1)*naxes(2)

   nbuf=yout

   j=0

   do while (npixels.gt.0)
      !read in 1 column at a time
      nbuffer=min(nbuf,npixels)

      j=j+1
      !write(0,*) j,xout,npixels
      !find max and min values
      do i=1,nbuffer
         !if ((j.eq.1595).and.(i.eq.194)) then
         !  write(0,*) j,i,k,pixels(j,i)/dngrpfac
         !endif
         buffer(i)=pixels(j,i)/dngrpfac !we are looping over naxis=2 to roate image.
      enddo
      call ftpprd(funit,group,firstpix,nbuffer,buffer,status)

      !update pointers and counters

      npixels=npixels-nbuffer
      firstpix=firstpix+nbuffer
   enddo

enddo


return
end subroutine writefitsdata
