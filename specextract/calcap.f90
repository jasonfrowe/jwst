program calcap
!Estimated aperture to enclose 90,95 and 99% of flux.
use precision
implicit none
integer :: iargc,nunit,filestatus,ntrace,nlines,nfit,i,j,nap,nkeysmax,  &
   nxmax,nymax,status,nkeys,dumi,i1,i2
integer, dimension(2) :: naxes
real(double) :: apstep,apmax,sq2pi,Pi,tflux,eflux,sqpid2,t1,t2,erf,r2,  &
   apinit,fpec,bpix,Rmin,Rmax,cen,tpflux
real(double), allocatable, dimension(:) :: psf,ap,aplevel,pap
real(double), allocatable, dimension(:,:) :: solpsf,Image,tImage
character(80) :: tracename,Imagename
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

!adjustable parameters
apstep=0.01 !stepsize to increase aperture.
apmax=100.0 !maximum aperture to accept
apinit=10.0  !initial aperture to try

!constants
Pi=acos(-1.d0)       !define Pi
sq2pi=sqrt(2.0d0*pi) !define sqrt(2*pi)
sqpid2=sqrt(Pi/2.0d0) !define sqrt(pi/2)
r2=sqrt(2.0d0)        !defined sqrt(2)

if(iargc().lt.2)then
   write(0,*) "Usage : calcap <Image> <tracepsf.dat>"
   write(0,*) "  <Image>        : FITS file containing SOSS data"
   write(0,*) "  <tracepsf.dat> : Output from spectextract"
   stop
endif

!get filename for trace
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

nap=3 !number of apertures to calculate
allocate(ap(nap),aplevel(nap),pap(nap))
aplevel(1)=0.90d0 !how much flux to encompass.
aplevel(2)=0.95d0
aplevel(3)=0.99d0

!234
!567
!890

do i=1,nlines
   !calculate total flux
   tflux=sq2pi*solpsf(i,8)*solpsf(i,2)*solpsf(i,4)+                     &
         sq2pi*solpsf(i,8)*solpsf(i,5)*solpsf(i,7)+                     &
         sq2pi*solpsf(i,8)*solpsf(i,10)

   do j=1,nap
      ap(j)=apinit
      fpec=0.0d0
      do while(fpec.le.aplevel(j))
         if(ap(j).gt.apmax) exit
         !calculate flux in aperture.
         t1=(solpsf(i,3)+ap(j))/(r2*solpsf(i,4))
         t2=(solpsf(i,3)-ap(j))/(r2*solpsf(i,4))
         eflux=sqpid2*solpsf(i,8)*solpsf(i,2)*solpsf(i,4)*(erf(t1)-erf(t2))
         t1=(solpsf(i,6)+ap(j))/(r2*solpsf(i,7))
         t2=(solpsf(i,6)-ap(j))/(r2*solpsf(i,7))
         eflux=eflux+sqpid2*solpsf(i,8)*solpsf(i,5)*solpsf(i,7)*(erf(t1)-erf(t2))
         t1=(+ap(j))/(r2*solpsf(i,10))
         t2=(-ap(j))/(r2*solpsf(i,10))
         eflux=eflux+sqpid2*solpsf(i,8)*solpsf(i,10)*(erf(t1)-erf(t2))

         fpec=eflux/tflux

!         write(0,*) "Tflux: ",i,eflux,eflux/tflux,ap(j)

         ap(j)=ap(j)+apstep
      enddo
      ap(j)=ap(j)-apstep
   enddo

   cen=solpsf(i,9)+(solpsf(i,3)+solpsf(i,6))/2.0d0
   tpflux=sum(Image(i,max(1,int(cen)-40):min(naxes(2),int(cen)+40)))

   do j=1,nap
      pap(j)=apinit
      fpec=0.0d0
      do while(fpec.le.aplevel(j))
         if(pap(j).gt.apmax) exit
         i1=max(1,int(cen-pap(j)))
         i2=min(naxes(2),int(cen+pap(j)))
!         write(0,*) i1,i2,cen
         eflux=sum(Image(i,i1:i2))

         fpec=eflux/tpflux

         pap(j)=pap(j)+1.0d0
      enddo

   enddo

   if(sum(ap).lt.apmax*3.0-1.0)then
      write(6,'(I4,6(1X,F6.2),1X,F11.2)') i,(2.0d0*ap(j),j=1,nap),      &
         (2.0d0*pap(j),j=1,nap),tflux
   endif

enddo


end program calcap
