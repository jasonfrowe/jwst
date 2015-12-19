program apextract
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: iargc,nkeysmax,nxmax,nymax,status,nkeys,dumi,filestatus,     &
   nunit,ntrace,nlines,nfit,i,j,k,nplot
integer, dimension(2) :: naxes
real, allocatable, dimension(:) :: px,py
real(double) :: bpix,Rmin,Rmax,tavg,avgsplitlength,tilt
real(double), allocatable, dimension(:) :: psf
real(double), allocatable, dimension(:,:) :: Image,tImage,solpsf,dTrace,&
   apfluxl,apfluxu
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
interface
   subroutine displayfits(nxmax,nymax,parray,bpix,tavg)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray
      real(double), intent(inout) :: bpix,tavg
   end subroutine displayfits
end interface
interface
   subroutine apsplit(naxes,Image,nlines,nTrace,solpsf,apfluxl,apfluxu)
      use precision
      implicit none
      integer, intent(in) :: nlines,nTrace
      integer, dimension(2), intent(in) :: naxes
      real(double), dimension(:,:), intent(inout) :: Image,solpsf,      &
         apfluxl,apfluxu
   end subroutine apsplit
end interface
interface
   subroutine xcorr(nlines,nTrace,apfluxl,apfluxu,avgsplitlength,tilt)
      use precision
      implicit none
      integer, intent(in) :: nlines,nTrace
      real(double), intent(in) :: avgsplitlength,tilt
      real(double), dimension(:,:), intent(in) :: apfluxl,apfluxu
   end subroutine xcorr
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

!Rotate image to deal with CV3
if(naxes(1).lt.naxes(2))then
   Image=transpose(Image)
   dumi=naxes(1)
   naxes(1)=naxes(2)
   naxes(2)=dumi
endif

write(0,*) "Image size: ",naxes

!resize image to minimize memory usage.
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

!open up pgplot window
call pgopen('/xserve')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,2)
call pgpage()

tavg=0.0 !displays a time on the image
call displayfits(nxmax,nymax,Image,bpix,tavg)

!read in trace
nunit=10 !unit number for file list
open(unit=nunit,file=tracename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",tracename
   stop
endif
read(nunit,*,iostat=filestatus) ntrace,nlines
if(filestatus.ne.0) then
   write(0,*) "Error reading line 1 : ",tracename
   stop
endif
if(nlines.ne.naxes(1))then
   write(0,*) "Image columns: ",naxes(1)
   write(0,*) "Trace columns: ",nlines
endif

nfit=1+9*ntrace
allocate(solpsf(nlines,nfit),psf(nfit))

do
   read(nunit,*,iostat=filestatus) i,(psf(j),j=1,nfit)
   if(filestatus == 0) then
      if(i.gt.nlines)then
         write(0,*) "Error: Trace-read has too many entries"
         write(0,*) "i     : ",i
         write(0,*) "nlines: ",nlines
         stop
      endif
      solpsf(i,:)=psf(:) !assign PSF model to array
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
write(0,*) "Trace PSF read ",i

!generate trace positions from PSF
allocate(dTrace(nlines,nTrace))
do i=1,nlines
   do k=1,ntrace
      dTrace(i,k)=solpsf(i,9+9*(k-1))+(solpsf(i,3+9*(k-1))+             &
         solpsf(i,6+9*(k-1)))/2.0d0
   enddo
enddo
!plot trace
allocate(px(naxes(1)),py(naxes(1)))
do i=1,nTrace
   nplot=0
   do j=1,naxes(1)
      if(dTrace(j,i).gt.0.1)then
         nplot=nplot+1
         px(nplot)=real(j)
         py(nplot)=real(dTrace(j,i))
      endif
   enddo
   call pgsci(2+i)
   call pgline(nplot,px,py)
enddo

allocate(apfluxl(nlines,nTrace),apfluxu(nlines,nTrace))
call apsplit(naxes,Image,nlines,nTrace,solpsf,apfluxl,apfluxu)

!estimate the average difference between the two peaks of the PSF
k=0
do i=1,nlines
   if( ( abs(solpsf(i,3)).gt.0.01 ).and.( abs(solpsf(i,3)).gt.0.01 ) ) then
      k=k+1
      avgsplitlength=avgsplitlength+solpsf(i,3)-solpsf(i,6)
   endif
enddo
avgsplitlength=avgsplitlength/dble(k)
!write(0,*) "avgsplitlength", avgsplitlength
call xcorr(nlines,nTrace,apfluxl,apfluxu,avgsplitlength,tilt)
write(6,'(A6,A32,1X,F6.3)') "Tilt: ",Imagename,tilt

call pgclos()

end program apextract
