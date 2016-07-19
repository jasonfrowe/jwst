program tiltest
!Jason Rowe 2016 - jasonfrowe@gmail.com
use precision
implicit none
integer :: iargc,nkeysmax,nxmax,nymax,status,nkeys,dumi,filestatus,     &
   nunit,ntrace,nlines,nfit,i,j,k,nplot,nfitline,nbin,nstart,nend,      &
   nlinesex,imax
integer, dimension(2) :: naxes
real, allocatable, dimension(:) :: px,py,px2,py2,bb
real(double) :: bpix,Rmin,Rmax,tavg,avgsplitlength,tilt
real(double), allocatable, dimension(:) :: psf
real(double), allocatable, dimension(:,:) :: Image,tImage,solpsf,dTrace,&
   apfluxl,apfluxu
character(80) :: Imagename,tracename,cline
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
   subroutine xcorr(nlines,nTrace,apfluxl,apfluxu,avgsplitlength,tilt,  &
    imax)
      use precision
      implicit none
      integer, intent(inout) :: nlines,nTrace,imax
      real(double), intent(inout) :: avgsplitlength,tilt
      real(double), dimension(:,:), intent(in) :: apfluxl,apfluxu
   end subroutine xcorr
end interface
interface
   subroutine fitline(nlines,nTrace,apfluxl,apfluxu,nfit)
      use precision
      implicit none
      integer, intent(in) :: nfit,nlines,nTrace
      real(double), dimension(:,:), intent(inout) :: apfluxl,apfluxu
   end subroutine fitline
end interface
interface
   subroutine polyfilter(nlines,nTrace,apfluxl,apfluxu,nbin)
      use precision
      implicit none
      integer, intent(in) :: nlines,nTrace,nbin
      real(double), dimension(:,:), intent(inout) :: apfluxl,apfluxu
   end subroutine polyfilter
end interface

nbin=20 !width of bin for filtering

!get filename
if(iargc().lt.2)then
   write(0,*) "Usage: tilt <Image> <tracepsf.dat> [nbin]"
   write(0,*) "  <Image>        : FITS file containing SOSS data"
   write(0,*) "  <tracepsf.dat> : Output from spectextract"
   write(0,*) "  [nbin]         : filter length (pixels, optional), default=20"
   write(0,*) "  [nstart]       : Column to start correlation with"
   write(0,*) "  [nend]         : Column to end correlation with"
   stop
endif
call getarg(1,Imagename)
call getarg(2,tracename)
if(iargc().ge.3)then
   call getarg(3,cline)
   read(cline,*) nbin
   if(nbin.le.1)then
      write(0,*) "nbin must be greater than 1"
      stop
   endif
endif

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

!now we grab commandline arguments and check against Image dimensions
nstart=1 !enables a portion of the spectrum to by analyzed - defaults
nend=2048
if(iargc().ge.4)then
   call getarg(4,cline)
   read(cline,*) nstart
   if(nstart.lt.1)then
      write(0,*) "nstart must be greater than 1"
      stop
   endif
endif
if(iargc().ge.5)then
   call getarg(5,cline)
   read(cline,*) nend
   if(nend.gt.naxes(1))then
      write(0,'(A26,I5)') "nend must be greater than ",naxes(1)
      stop
   endif
   if(nend.le.nstart)then
      write(0,*) "nend must be greater than nstart"
      stop
   endif
endif
nlinesex=nend-nstart+1 !number of points to extract.
if(nlinesex.lt.3)then
   write(0,*) "Must hvae nend - mstart > 3"
   stop
endif

!open up pgplot window
!call pgopen('/xserve')
call pgopen('?')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage()

tavg=0.0 !displays a time on the image
call displayfits(nxmax,nymax,Image,bpix,tavg,0.0d0)

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
!read in data from file line by line
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
deallocate(px,py)



!extract flux on each side of the trace via PSF model
allocate(apfluxl(nlines,nTrace),apfluxu(nlines,nTrace))
call apsplit(naxes,Image,nlines,nTrace,solpsf,apfluxl,apfluxu)

!fitting a polynomial to remove the continuum
!nfitline=12
!call fitline(nlines,nTrace,apfluxl,apfluxu,nfitline)

!using a polynomial based bandpass filter to remove the continuum
!nbin controls the filter
if(nbin.lt.nlines)then
   call polyfilter(nlines,nTrace,apfluxl,apfluxu,nbin)
endif

allocate(px(nlines),py(nlines),bb(4))
nplot=0
k=1 !plot first Trace
do i=1,nlines
!   if(apfluxl(i,k).gt.0.0)then
      nplot=nplot+1
      px(nplot)=real(i)
      py(nplot)=real(apfluxl(i,k))
!      write(0,*) px(nplot),py(nplot)
!   endif
enddo
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
bb(1)=minval(px(1:nplot))
bb(2)=maxval(px(1:nplot))
bb(3)=minval(py(1:nplot))
bb(4)=maxval(py(1:nplot))
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
!call pgwindow(50.0,250.0,minval(py),maxval(py(50:250)))

allocate(px2(4),py2(4))
px2(1)=nstart
py2(1)=minval(py(1:nplot))
px2(2)=nstart
py2(2)=maxval(py(1:nplot))
px2(3)=nend
py2(3)=py2(2)
px2(4)=nend
py2(4)=py2(1)
call pgsci(4)
call pgsfs(1)
call pgpoly(4,px2,py2)
call pgsci(1)
deallocate(px2,py2)

call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "X (pixels)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Y (counts)")
call pgline(nplot,px,py)

nplot=0
k=1 !plot first Trace
do i=1,nlines
!   if(apfluxu(i,k).gt.0.0)then
      nplot=nplot+1
      px(nplot)=real(i)
      py(nplot)=real(apfluxu(i,k))
!   endif
enddo
call pgpage()
call pgsci(1)
call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
bb(1)=minval(px(1:nplot))
bb(2)=maxval(px(1:nplot))
bb(3)=minval(py(1:nplot))
bb(4)=maxval(py(1:nplot))
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
!call pgwindow(50.0,250.0,minval(py),maxval(py(50:250)))

allocate(px2(4),py2(4))
px2(1)=nstart
py2(1)=minval(py(1:nplot))
px2(2)=nstart
py2(2)=maxval(py(1:nplot))
px2(3)=nend
py2(3)=py2(2)
px2(4)=nend
py2(4)=py2(1)
call pgsci(4)
call pgsfs(1)
call pgpoly(4,px2,py2)
call pgsci(1)
deallocate(px2,py2)

call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.16*(bb(4)-bb(3)),0.0,0.5,         &
   "X (pixels)")
call pgptxt(bb(1)-0.04*(bb(2)-bb(1)),(bb(4)+bb(3))/2,90.0,0.5,          &
   "Y (counts)")
call pgline(nplot,px,py)
deallocate(px,py)


!estimate the average difference between the two peaks of the PSF
!this length is used to calculate the angle.
k=0
do i=1,nlines
   if( ( abs(solpsf(i,3)).gt.0.01 ).and.( abs(solpsf(i,3)).gt.0.01 ) ) then
      k=k+1
      avgsplitlength=avgsplitlength+solpsf(i,3)-solpsf(i,6)
   endif
enddo
avgsplitlength=avgsplitlength/dble(k)
write(0,*) "avgsplitlength", avgsplitlength

!calculate cross-correlation
!call xcorr(nlines,nTrace,apfluxl,apfluxu,avgsplitlength,tilt,imax)
call xcorr(nlinesex,nTrace,apfluxl(nstart:nend,:),                      &
   apfluxu(nstart:nend,:),avgsplitlength,tilt,imax)
write(6,'(A6,A32,1X,F9.4,1X,I5)') "Tilt: ",Imagename,tilt,imax

call pgclos()

end program tiltest
