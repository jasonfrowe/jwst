subroutine readcfits(filename,naxes,imax,Imagecube,bpix)
use precision
implicit none
integer, dimension(:) :: naxes,imax
real(double) :: bpix
real(double), dimension(:,:,:) :: Imagecube
character(200) :: filename
!local vars
integer :: istatus,unitfits,readwrite,dumi,nkeys,nspace,i,nfound,group, &
 firstpix,nbuf,j,naxis
real(double) :: nullval
real(double), allocatable, dimension(:) :: buffer
character(80) :: record
character(80), allocatable, dimension(:) :: header
logical :: anynull

! status will report errors.  No errors means status=0.
! initalize value of status
istatus=0
! gets an unused unit number to open fits file
call ftgiou(unitfits,istatus)
! setting to zero makes fits file readwrite
readwrite=0
! open this fits file
call ftopen(unitfits,filename,readwrite,dumi,istatus)
if(istatus.ne.0)then
   write(0,*) "Status: ",istatus
   write(0,*) "Cannot open "
   write(0,'(A80)') filename
   stop
endif

nkeys=0
! get number of headers in image
call ftghsp(unitfits,nkeys,nspace,istatus)
allocate(header(nkeys))
do i=1,nkeys
   call ftgrec(unitfits,i,record,istatus)
   header(i)=record
!   write(6,'(A80)') header(i)
enddo

!read in dimension of image
call ftgidm(unitfits,naxis,istatus)
!write(0,*) "naxis: ",naxis

call ftgknj(unitfits,'NAXIS',1,3,naxes,nfound,istatus)
if((naxes(1).gt.size(Imagecube(:,1,1))).or.                             &
 (naxes(2).gt.size(Imagecube(1,:,1))).or.                               &
 (naxes(3).gt.size(Imagecube(1,1,:))))then
   write(0,*) "inadequate space for FITS."
   write(0,*) "Needed: ", naxes(1),naxes(2),naxes(3)
   write(0,*) "Available: ", size(Imagecube(:,1,1)),                    &
      size(Imagecube(1,:,1)),size(Imagecube(1,1,:))
   stop
endif

!Check that it found both NAXIS1 and NAXIS2 keywords.
if (nfound.ne.3)then
   write(6,*) 'READIMAGE failed to read the NAXISn keywords.'
   stop
endif

write(0,*) "Dimensions ",naxes(1),naxes(2),naxes(3)

!npixels=naxes(1)*naxes(2)*naxes(3)
group=1
firstpix=1
nullval=bpix
nbuf=naxes(1)
allocate(buffer(nbuf))


do j=1,naxes(3)
   do i=1,naxes(2)
      call ftgpvd(unitfits,group,firstpix,nbuf,nullval,buffer,          &
       anynull,istatus)
      Imagecube(1:nbuf,i,j)=buffer(1:nbuf)
      if(istatus.ne.0)then
         write(0,*) "istatus: ",istatus
         write(0,*) "Error reading in Image"
      endif
      firstpix=firstpix+nbuf
   enddo
enddo

!write(0,*) "min/max: ",minval(Imagecube(1:naxes(1),1:naxes(2),1:naxes(3))),&
!   maxval(Imagecube(1:naxes(1),1:naxes(2),1:naxes(3)))

!do i=1,naxes(3)
!   write(6,*) Imagecube(200,100,i)
!   write(0,*) sum(Imagecube(1:naxes(1),1:naxes(2),i))
!enddo

return
end subroutine readcfits
