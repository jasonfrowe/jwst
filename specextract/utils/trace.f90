subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf)
!generates trace
use precision
implicit none
integer :: nkeys,nkeysmax,nKsize,i,nline,nksd2,xm,xp,j,nmaxf,k,nTrace,  &
   nTp,jj,ncut,ncutpsf,ncutd2
integer, dimension(2) :: naxes,knaxes
integer, allocatable, dimension(:) :: isol,isolfirst
real(double) :: Kmin,Kmax,bpix,maxf,S,Sxy,Sxx,Sx,Sy,df,bcut,f,          &
   triplegaussian,sq2pi,pi,tracetest,tcut
real(double), dimension(:,:) :: Image,dTrace,bf
real(double), allocatable, dimension(:) :: line,lpsf,lpsftemp,psfwork,  &
   sol,model,solnew,amp,solfirst
real(double), allocatable, dimension(:,:) :: Kernel,psf
character(80) :: kfile
character(80), allocatable, dimension(:) :: header
!plotting
real :: pmin
real, allocatable, dimension(:) :: px,py

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
   subroutine loadPSFinit(ntrace,sol,ncutpsf,nline,dtrace,line)
      use precision
      implicit none
      integer, intent(inout) :: ntrace,ncutpsf,nline
      real(double), dimension(:), intent(inout) :: sol,line
      real(double), dimension(:,:), intent(inout) :: dtrace
   end subroutine loadPSFinit
end interface
interface
   subroutine psfmodel1d(npt,model,ntrace,sol)
      use precision
      implicit none
      integer, intent(in) :: npt,ntrace
      real(double), dimension(:), intent(in) :: sol
      real(double), dimension(:), intent(inout) :: model
   end subroutine psfmodel1d
end interface
interface
   subroutine modelline(npt,line,ntrace,sol,isol)
      use precision
      implicit none
      integer, intent(in) :: npt,ntrace
      integer, dimension(:), intent(inout) :: isol
      real(double), dimension(:), intent(inout) :: line,sol
   end subroutine modelline
end interface


!constants
nKsize=64        !size of Kernel for import (note.. this gets resized)
nksd2=nKsize/2 !precompute size of PSF/2
Pi=acos(-1.d0)       !define Pi
sq2pi=sqrt(2.0d0*pi) !define sqrt(2*pi)
!parameters to control trace
ncut=35          !width of spectrum to zero out
bcut=10.0d0      !threshold for finding a trace
tcut=10.0d0       !threshold for traces to jump
ncutpsf=24       !width of PSF to fit - must be less than nKsize/2
if(ncutpsf.gt.nKsize)then  !check ncutpsf value is valid
   write(0,*) "Error: ncutpsf must be less than nKsize"
   write(0,*) "ncutpsf : ",ncutpsf
   write(0,*) "nKsize: ",nKsize
   stop !if ncutpsf is invalid then stop prog.
endif
ncutd2=ncutpsf/2 !precompute ncutd2

!read in a Kernel to be used as a guess for the PSF
kfile="Kernels1/psf_2100nm_x10_oversampled.fits"
nkeysmax=800 !maximum number of lines in the header we will accept.
allocate(header(nkeysmax)) !allocate space for header (we really don't need it)
allocate(Kernel(nKsize,nKsize))  !allocate space for Kernel
call getfits(kfile,knaxes,Kernel,Kmin,Kmax,nkeys,header,bpix) !read in FITS
write(0,*) "knaxes: ",knaxes

!generate 1-D PSF for cross-correlation
allocate(lpsf(nKsize),psfwork(nKsize))
lpsf=Sum(Kernel,2)
lpsf=lpsf(nKsize:1:-1) !looks like the PSF needs to be flipped to match spgen
!lpsf=lpsf-minval(lpsf) !remove zero point?

allocate(psf(nKsize,ntrace))
do i=1,nTrace    !make a trace for each order.
   psf(:,i)=lpsf(:)  !our initial guess for each order to be refined
enddo

!allocate space for line-scan
allocate(line(naxes(2)))  !this will contain each column of data to scan

line=Image(nline,:)    !we start at user-defined 'nline'
!subtract off minimum value - maybe a true sky value would be better
line=line-minval(line)

do k=1,nTrace !loop over expected number of traces

   maxf=0.0d0 !initalize scan for maximum value of f
   do i=1,naxes(2) !scan over the line
   !xm, xp give position along line to correlate with PSF
      xm=max(1,i-ncutd2+1)
      xp=min(naxes(2),i+ncutd2)
      !calculate the correlation, f = line * PSF centered at desired position
      f=Sum(line(xm:xp)*lpsf(nksd2-(i-xm):nksd2+(xp-i)))

      !when f is maximized we have a best-match to the PSF
      if(f.gt.maxf)then
         maxf=f
         nmaxf=i   !store position of maximum 'f' value
      endif
   enddo
   write(0,*) "maxf: ",nmaxf,maxf
   dTrace(nline,k)=nmaxf

   !position of trace
   i=nmaxf
   xm=max(1,i-ncutd2+1)      !we will now fit the PSF to the line
   xp=min(naxes(2),i+ncutd2)
   !calculate best-fit scaling of PSF to line.
   S=1.0d0
   Sxy=Sum(line(xm:xp)*lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sxx=Sum(lpsf(nksd2-(i-xm):nksd2+(xp-i))*lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sx =Sum(lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sy =Sum(line(xm:xp))
   bf(nline,k)=(S*Sxy-Sx*Sy)/(S*Sxx-Sx*Sx)  ! line = psf * bf
   write(0,*) "bf: ",bf(nline,k)
!  calculate residual
   line(xm:xp)=line(xm:xp)-lpsf(nksd2-(i-xm):nksd2+(xp-i))*bf(nline,k)

!  plot the fitted PSF against the line
!  PGPLOT uses REAL*4, so we convert dble -> real
   allocate(px(size(line(xm:xp))),py(size(line(xm:xp))))
   do j=xm,xp
      px(j-xm+1)=real(j)  !X axis
   enddo
   py=real(lpsf(nksd2-(i-xm):nksd2+(xp-i))*bf(nline,k)) !Y-axis
!   pmin=minval(py)
!   py=log10(py-pmin+1.0d0)
   call pgsci(2+k)  !change plotting colour
   call pgline(size(line(xm:xp)),px,py) !plot a line
   call pgsci(1) !change plotting colour back to default
   deallocate(px,py) !de-allocate plotting variables

!There can be siginficant differences compare to first pass, so
!zero out residuals to avoid problems.
   xm=max(1,i-ncut+1)
   xp=min(naxes(2),i+ncut)
   do j=xm,xp
      line(j)=0.0d0 !zero out value - should actually be set to sky.
   enddo
enddo


!Lets model the line with a PSF.
line=Image(nline,:)    !we start at user-defined 'nline'
!load initial-Guess in sol.  The parameter isol controls with variables
!are fit
allocate(sol(ntrace*9+1),isol(ntrace*9+1),solnew(ntrace*9+1),amp(ntrace))
call loadPSFinit(ntrace,sol,ncutpsf,nline,dtrace,line)
!write(0,*) "****",(sol(i),i=1,ntrace*9+1)
write(0,'(A2,1X,I4,3(1X,F11.3))') "*1",nline,                           &
   (sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0,k=1,3)
isol=1  !fit all variables
isol(1)=0 !do not fit zero line
call modelline(naxes(2),line,ntrace,sol,isol)
do k=1,ntrace
   amp(k)=sq2pi*solnew(8+9*(k-1))*solnew(2+9*(k-1))*solnew(4+9*(k-1))+  &
          sq2pi*solnew(8+9*(k-1))*solnew(5+9*(k-1))*solnew(7+9*(k-1))+  &
          sq2pi*solnew(8+9*(k-1))*solnew(10+9*(k-1))
   bf(i,k)=amp(k)
enddo
do k=1,ntrace
!   dTrace(nline,k)=sol(9+9*(k-1))
   dTrace(nline,k)=sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0
enddo
write(0,'(A2,1X,I4,3(1X,F11.3))') "**",nline,                           &
   (sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0,k=1,3)

allocate(model(size(line)))
allocate(px(size(line)),py(size(line)))
do j=1,size(line)
   px(j)=real(j)  !X axis
enddo
call psfmodel1d(size(line),model,ntrace,sol)
py=real(model) !Y-axis
!   pmin=minval(py)
!   py=log10(py-pmin+1.0d0)
call pgsci(6)  !change plotting colour
call pgline(size(line),px,py) !plot a line
call pgsci(1) !change plotting colour back to default
deallocate(px,py) !de-allocate plotting variables
deallocate(model)
!write(0,*) "pause.."
!read(5,*)

!fixing the shape of the PSF model (only central position and amplitude)
do k=1,ntrace
   isol(2+9*(k-1))=0
   isol(3+9*(k-1))=0
   isol(4+9*(k-1))=0
   isol(5+9*(k-1))=0
   isol(6+9*(k-1))=0
   isol(7+9*(k-1))=0
   isol(10+9*(k-1))=0
enddo

allocate(solfirst(ntrace*9+1),isolfirst(ntrace*9+1))
solfirst=sol !save solution from first line
isolfirst=isol

!Next part is to use the found positions and develop the trace.

do i=nline+1,naxes(1) !forward direction
   line=Image(i,:)  !get next column from image to use.
   line=line-minval(line) !offset to zero - change to sky-subtraction
   solnew=sol !save current solution
   call modelline(naxes(2),line,ntrace,solnew,isol) !model
!  now check the amplitudes and changes in trace
   do k=1,ntrace
      amp(k)=sq2pi*solnew(8+9*(k-1))*solnew(2+9*(k-1))*solnew(4+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(5+9*(k-1))*solnew(7+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(10+9*(k-1))
      bf(i,k)=amp(k)
      if(amp(k).lt.bcut)then  !if amplitude is too low - kill trace
         do j=2+9*(k-1),10+9*(k-1)
            solnew(j)=0.0d0
            isol(j)=0
            bf(i,k)=0.0d0
         enddo
      endif
      tracetest=sol(9+9*(k-1))+(solnew(3+9*(k-1))+solnew(6+9*(k-1)))/2.0d0
      if(abs(dTrace(i-1,k)-tracetest).gt.tcut)then !if trace jumps. kill trace
         write(0,*) "Trace Jumped ",k,tracetest,dTrace(i-1,k)
         do j=2+9*(k-1),10+9*(k-1)
            solnew(j)=0.0d0
            isol(j)=0
            bf(i,k)=0.0d0
         enddo
      endif
   enddo
!   write(0,*) "amp: ",(amp(k),k=1,3)
   sol=solnew
   do k=1,ntrace
!      dTrace(i,k)=sol(9+9*(k-1))
      dTrace(i,k)=sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0
   enddo

   write(0,'(A2,1X,I4,3(1X,F11.3))') "**",i,(dTrace(i,k),k=1,3)

enddo

!restore solution from nline.
sol=solfirst
isol=isolfirst
do i=nline-1,1,-1 !negative direction
   line=Image(i,:)
   line=line-minval(line)
   solnew=sol !save current solution
   call modelline(naxes(2),line,ntrace,solnew,isol) !model
!  now check the amplitudes and changes in trace
   do k=1,ntrace
      amp(k)=sq2pi*solnew(8+9*(k-1))*solnew(2+9*(k-1))*solnew(4+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(5+9*(k-1))*solnew(7+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(10+9*(k-1))
      bf(i,k)=amp(k)
      if(amp(k).lt.bcut)then  !amplitude is too low.. kill trace
         do j=2+9*(k-1),10+9*(k-1)
            solnew(j)=0.0d0   !set PSF model parameter to zero
            isol(j)=0         !disable variables
            bf(i,k)=0.0d0
         enddo
      endif
      tracetest=sol(9+9*(k-1))+(solnew(3+9*(k-1))+solnew(6+9*(k-1)))/2.0d0
      if(abs(dTrace(i+1,k)-tracetest).gt.tcut)then !trace jumped alot
         write(0,*) "Trace Jumped ",k,tracetest,dTrace(i+1,k)
         do j=2+9*(k-1),10+9*(k-1)  !kill trace
            solnew(j)=0.0d0
            isol(j)=0
            bf(i,k)=0.0d0
         enddo
      endif
   enddo
!   write(0,*) "amp: ",(amp(k),k=1,3)
   sol=solnew
   do k=1,ntrace
!      dTrace(i,k)=sol(9+9*(k-1))
      dTrace(i,k)=sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0
   enddo

   write(0,'(A2,1X,I4,3(1X,F11.3))') "**",i,(dTrace(i,k),k=1,3)

enddo


write(0,*) "End of Trace.. "
return
end subroutine trace
