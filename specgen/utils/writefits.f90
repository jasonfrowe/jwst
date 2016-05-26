subroutine writefits(nxmax,nymax,parray,fileout,time)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nxmax,nymax,nkeys,nstep,status,blocksize,bitpix,naxis,funit, &
   npixels,group,firstpix,nbuf,i,j,nbuffer
integer, dimension(2) :: naxes
integer, dimension(4) :: nr
real(double) :: time
real(double), allocatable, dimension(:) :: buffer
real(double), dimension(:,:) :: parray
character(80) :: fileout,record
logical simple,extend

naxes(1)=nxmax !size of image to write to FITS file
naxes(2)=nymax

status=0
!if file already exists.. delete it.
call deletefile(fileout,status)
!get a unit number
call ftgiou(funit,status)
!Create the new empty FITS file.  The blocksize parameter is a
!historical artifact and the value is ignored by FITSIO.
blocksize=1
status=0
call ftinit(funit,fileout,blocksize,status)
if(status.ne.0)then
   write(0,*) "Status: ",status
   write(0,*) "Critial Error open FITS for writing"
   write(0,'(A80)') fileout
endif

!Initialize parameters about the FITS image.
!BITPIX = 16 means that the image pixels will consist of 16-bit
!integers.  The size of the image is given by the NAXES values.
!The EXTEND = TRUE parameter indicates that the FITS file
!may contain extensions following the primary array.
simple=.true.
bitpix=-32
naxis=2
extend=.true.

!Write the required header keywords to the file
call ftphpr(funit,simple,bitpix,naxis,naxes,0,1,extend,status)

!Write the array to the FITS file.
npixels=naxes(1)*naxes(2)
group=1
firstpix=1
nbuf=naxes(1)
j=0

allocate(buffer(nbuf))
do while (npixels.gt.0)
!read in 1 column at a time
   nbuffer=min(nbuf,npixels)

   j=j+1
!find max and min values
   do i=1,nbuffer
      buffer(i)=parray(i,j)
   enddo

   call ftpprd(funit,group,firstpix,nbuffer,buffer,status)

!update pointers and counters

   npixels=npixels-nbuffer
   firstpix=firstpix+nbuffer

enddo

!write(6,*) "ftprec:",status
write(record,'(A8,A3,F12.8)') 'HJD     ','=  ',time
write(6,'(a80)') record
call ftprec(funit,record,status)
!!write(record,'(A8,A3,F10.1)') 'DATAMIN ','=  ',-10000.0
!write(6,'(a80)') record
!!call ftprec(funit,record,status)
!write(6,*) "ftprec:",status


!close fits file
call ftclos(funit,status)
call ftfiou(funit,status)

return
end subroutine writefits
