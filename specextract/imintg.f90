program imintg
!reads in data cubes and integrates.
use precision
implicit none
integer :: iargc,i,j,npt,k,nplot,ipflag
integer, allocatable, dimension(:) :: naxes,imax
real :: t1,t2
real, allocatable, dimension(:) :: bb,px,py,lx,ly
real(double) :: bpix,sat,a,b,abdev,dnaxes3
real(double), allocatable, dimension(:) :: x,y
real(double), allocatable, dimension(:,:) :: Image
real(double), allocatable, dimension(:,:,:) :: Imagecube,tImagecube
character(200) :: filename,fileout

interface
   subroutine readcfits(filename,naxes,imax,Imagecube,bpix)
      use precision
      implicit none
      integer, dimension(:) :: naxes,imax
      real(double) :: bpix
      real(double), dimension(:,:,:) :: Imagecube
      character(200) :: filename
   end subroutine
end interface
interface
   subroutine writefits(nxmax,nymax,parray,fileout)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray
      character(80) :: fileout
   end subroutine
end interface

!parameters controlling algorithms
sat=65535.0d0  !saturation value
bpix=-1.0d10 !marking bad pixels
ipflag=1 !0=plot, 1=no plot

if(iargc().lt.2)then
   write(0,*) "Usage: imintg Datacube Image"
   stop
endif

call getarg(1,filename)
call getarg(2,fileout)

allocate(naxes(3),imax(3))
imax(1)=2048
imax(2)=2048
imax(3)=100
allocate(Imagecube(imax(1),imax(2),imax(3)))
call readcfits(filename,naxes,imax,Imagecube,bpix)

!compact memory requirements for array
if((naxes(1).ne.imax(1)).or.(naxes(2).ne.imax(2)).or.(naxes(3).ne.imax(3)))then
   allocate(tImagecube(naxes(1),naxes(2),naxes(3)))
   tImagecube(1:naxes(1),1:naxes(2),1:naxes(3))=Imagecube(1:naxes(1),1:naxes(2),1:naxes(3))
   deallocate(Imagecube)
   imax(1)=naxes(1)
   imax(2)=naxes(2)
   imax(3)=naxes(3)
   allocate(Imagecube(imax(1),imax(2),imax(3)))
   Imagecube=tImagecube
   deallocate(tImagecube)
endif

if(ipflag.eq.0)then
   call pgopen('?')
   call PGPAP (8.0 ,1.0) !use a square 8" across
   !call pgsubp(1,4)
   call pgpage()
   call pgslw(3) !thicker lines
   !call pgask(.false.) !whether to prompt when changing frames
endif

allocate(x(imax(3)),y(imax(3)),Image(naxes(1),naxes(2)))
if(ipflag.eq.0) allocate(bb(4),px(imax(3)),py(imax(3)),lx(2),ly(2))

dnaxes3=dble(naxes(3))

!loop over 2D image to measure ramps
do i=1,naxes(1)
   do j=1,naxes(2)
      npt=0 !initialize number of points to fit to zero.
      do k=1,naxes(3) !collect all values for a single pixel
         if((Imagecube(i,j,k).lt.sat).and.(Imagecube(i,j,k).gt.bpix))then
            npt=npt+1
            x(npt)=dble(npt)
            y(npt)=Imagecube(i,j,k)
         endif
      enddo

!     robust fitting
      if(npt.ge.2)then
         call medfit(x,y,npt,a,b,abdev)
!        create Image
         Image(i,j)=b*dnaxes3 !construct image
      else
         Image(i,j)=Imagecube(i,j,naxes(3)) !construct image
      endif

!     Plotting
      if(ipflag.eq.0)then
         nplot=npt
         px(1:npt)=real(x(1:npt))
         py(1:npt)=real(y(1:npt))
         bb(1)=0.0
         bb(2)=maxval(px(1:nplot))+1.0
         t1=minval(py(1:nplot))
         t2=maxval(py(1:nplot))
         bb(3)=t1-(t2-t1)*0.10
         bb(4)=t2+(t2-t1)*0.10
         call pgwindow(bb(1),bb(2),bb(3),bb(4))
         CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
         call pgpt(nplot,px,py,17)
!        plotting solution
         lx(1)=bb(1)
         lx(2)=bb(2)
         ly(1)=real(a)+real(b)*lx(1)
         ly(2)=real(a)+real(b)*lx(2)
         call pgsci(2)
         call pgline(2,lx,ly)
         call pgsci(1)
         call pgpage()
      endif

   enddo
enddo

if(ipflag.eq.0) call pgclos()

call writefits(naxes(1),naxes(2),Image,fileout)

end program imintg
