subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf,solpsf,posguess)
!Jason Rowe 2015 - jasonfrowe@gmail.com
!generates trace - Version 2.0
use precision
implicit none
integer, parameter :: nKsize=64
integer :: nkeys,nkeysmax,i,nline,nksd2,xm,xp,j,nmaxf,k,nTrace,  &
   nTp,jj,ncut,ncutpsf,ncutd2,ilinkpsf,nbin
integer, dimension(2) :: naxes,knaxes
integer, allocatable, dimension(:) :: isol,isolfirst
real(double) :: Kmin,Kmax,bpix,maxf,S,Sxy,Sxx,Sx,Sy,df,bcut,f,          &
   triplegaussian,sq2pi,pi,tracetest,tcut,sky,std,snr,signal
real(double), dimension(:) :: posguess
real(double), dimension(:,:) :: Image,dTrace,bf,solpsf
real(double), dimension(nKsize) :: lpsf
real(double), allocatable, dimension(:) :: line,lpsftemp,psfwork,  &
   sol,model,solnew,amp,solfirst,solbin
real(double), allocatable, dimension(:,:) :: Kernel,psf
character(80) :: kfile
character(80), allocatable, dimension(:) :: header
!plotting
real :: pmin
real, allocatable, dimension(:) :: px,py
data  lpsf /3.746300735497921E-006,7.640107659745432E-006,4.542999963486538E-006, &
  4.929265102315838E-006,  9.714837387708730E-006,  5.671904909021475E-006, &
  4.548023730510664E-006,  1.022713226439542E-005,  6.886893882507295E-006, &
  8.177790144225927E-006,  1.357109057534278E-005,  8.916710340478584E-006, &
  1.239566539967818E-005,  2.781745489985332E-005,  2.509716449416999E-005, &
  2.631011432652208E-005,  4.830151269574756E-005,  5.380108778450451E-005, &
  7.547263667725956E-005,  1.022118883162726E-004,  1.420077972523748E-004, &
  2.362241206079752E-004,  3.385566380821325E-004,  5.477043893594158E-004, &
  6.696818098559376E-004,  5.493319867611035E-004,  4.720754680409556E-004, &
  2.991750642213908E-004,  3.058475983204190E-004,  3.109660592775787E-004, &
  2.226914950899106E-004,  2.979418360802288E-004,  3.397708704659941E-004, &
  2.990017218531538E-004,  2.758223087866440E-004,  3.294162992516503E-004, &
  2.381257536346881E-004,  3.407167609725814E-004,  5.361983993812380E-004, &
  6.230353641937803E-004,  6.140843798414508E-004,  5.070604643273580E-004, &
  3.306009460586345E-004,  2.371751859966409E-004,  1.155928608405077E-004, &
  6.370671124544813E-005,  7.242587988226523E-005,  3.417951946560471E-005, &
  2.799611461752616E-005,  3.403609616187131E-005,  1.659919620922157E-005, &
  1.873137450902895E-005,  1.825263581423098E-005,  9.144579730557822E-006, &
  8.962197003747896E-006,  9.059127708432868E-006,  5.408733372069818E-006, &
  7.623564329893584E-006,  7.990429056435600E-006,  4.048852234816991E-006, &
  6.761537376720472E-006,  7.205561307621622E-006,  4.591548022681025E-006, &
  4.112999194311184E-006/

interface
!   subroutine getfits(Refname,naxes,Ref,Rmin,Rmax,nkeys,header,bpix)
!      use precision
!      implicit none
!      integer :: nkeys
!      integer, dimension(2), intent(inout) :: naxes
!      real(double), intent(inout) :: Rmin,Rmax,bpix
!      real(double), dimension(:,:), intent(inout) :: Ref
!      character(80), intent(inout) :: Refname
!      character(80), dimension(:), intent(inout) :: header
!   end subroutine getfits
   subroutine loadPSFinit(ntrace,sol,ncutpsf,nline,dtrace,line)
      use precision
      implicit none
      integer, intent(inout) :: ntrace,ncutpsf,nline
      real(double), dimension(:), intent(inout) :: sol,line
      real(double), dimension(:,:), intent(inout) :: dtrace
   end subroutine loadPSFinit
   subroutine psfmodel1d(npt,model,ntrace,sol)
      use precision
      implicit none
      integer, intent(in) :: npt,ntrace
      real(double), dimension(:), intent(in) :: sol
      real(double), dimension(:), intent(inout) :: model
   end subroutine psfmodel1d
   subroutine modelline(npt,line,ntrace,sol,isol,ilinkpsf,nline)
      use precision
      implicit none
      integer, intent(inout) :: npt,ntrace,ilinkpsf,nline
      integer, dimension(:), intent(inout) :: isol
      real(double), dimension(:), intent(inout) :: line,sol
   end subroutine modelline
   function getsky(naxes,Image,std)
      use precision
      implicit none
      integer, dimension(2), intent(in) :: naxes
      real(double) :: getsky,std
      real(double), dimension(:,:), intent(in) :: Image
   end function getsky
end interface

!constants
!nKsize=64        !size of Kernel for import (note.. this gets resized)
nksd2=nKsize/2 !precompute size of PSF/2
Pi=acos(-1.d0)       !define Pi
sq2pi=sqrt(2.0d0*pi) !define sqrt(2*pi)
!parameters to control trace
ncut=35          !width of spectrum to zero out
bcut=3.0d0       !S/N threshold for finding a trace
tcut=10.0d0      !threshold for traces to jump
ncutpsf=24       !width of initial PSF to fit - must be less than nKsize/2
if(ncutpsf.gt.nKsize)then  !check ncutpsf value is valid
   write(0,*) "Error: ncutpsf must be less than nKsize"
   write(0,*) "ncutpsf : ",ncutpsf
   write(0,*) "nKsize: ",nKsize
   stop !if ncutpsf is invalid then stop prog.
endif
ncutd2=ncutpsf/2 !precompute ncutd2

!!read in a Kernel to be used as a guess for the PSF
!kfile="Kernels1/psf_2100nm_x10_oversampled.fits"
!nkeysmax=800 !maximum number of lines in the header we will accept.
!allocate(header(nkeysmax)) !allocate space for header (we really don't need it)
!allocate(Kernel(nKsize,nKsize))  !allocate space for Kernel
!call getfits(kfile,knaxes,Kernel,Kmin,Kmax,nkeys,header,bpix) !read in FITS
!write(0,*) "knaxes: ",knaxes
!
!!generate 1-D PSF for cross-correlation
!allocate(lpsf(nKsize),psfwork(nKsize))
!lpsf=Sum(Kernel,2)
!lpsf=lpsf(nKsize:1:-1) !looks like the PSF needs to be flipped to match spgen

allocate(psf(nKsize,ntrace))
do i=1,nTrace    !make a trace for each order.
   psf(:,i)=lpsf(:)  !our initial guess for each order to be refined
enddo

!allocate space for line-scan
allocate(line(naxes(2)))  !this will contain each column of data to scan

line=Image(nline,:)    !we start at user-defined 'nline'
sky=getsky(naxes,Image,std) !return sky and standard-deviation
write(0,*) "Sky: ",sky,std
!subtract off minimum value - maybe a true sky value would be better
line=line-sky


!this part is to get the initial positions of the traces for each order
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
!   call pgsci(2+k)  !change plotting colour
!   call pgline(size(line(xm:xp)),px,py) !plot a line
!   call pgsci(1) !change plotting colour back to default
   deallocate(px,py) !de-allocate plotting variables

!  calculate S/N
   signal=Sum(lpsf(nksd2-(i-xm):nksd2+(xp-i))*bf(nline,k))
   snr=signal/(sqrt(dble(xp-xm+1)*std))
   write(0,*) "SNR: ",snr

!  copy S/N into bf array to return
   bf(nline,k)=snr

!There can be siginficant differences compare to first pass, so
!zero out residuals to avoid problems.
   xm=max(1,i-ncut+1)
   xp=min(naxes(2),i+ncut)
   do j=xm,xp
      line(j)=0.0d0 !zero out value - should actually be set to sky.
   enddo
enddo

!now we have initial positions, if you want to override, then set
!dtrace(nline,k) to appropriate guess of position (k=trace#)
!e.g., dtrace(nline,1)=194
!We also can use commandline parameters for this:
do i=1,nTrace
   if(posguess(i).gt.0.0d0)then
      dtrace(nline,i)=posguess(i)
   endif
enddo


!Lets model the line with a PSF.
line=Image(nline,:)    !we start at user-defined 'nline'
line=line-sky
!load initial-Guess in sol.  The parameter isol controls with variables
!are fit
allocate(sol(ntrace*9+1),isol(ntrace*9+1),solnew(ntrace*9+1),amp(ntrace))
allocate(solbin(ntrace*9+1))
call loadPSFinit(ntrace,sol,ncutpsf,nline,dtrace,line)
!write(0,'(A2,1X,I4,3(1X,F11.3))') "*1",nline,                           &
!   (sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0,k=1,ntrace)
isol=1  !fit all variables
isol(1)=0 !do not fit zero line
ilinkpsf=0 !(0) - each order has own PSF, (1) PSF shapes linked
!This is the FIRST psf model update.. so keep the parameters free as
!possible.
call modelline(naxes(2),line,ntrace,sol,isol,ilinkpsf,nline)
do k=1,ntrace
   amp(k)=sq2pi*sol(8+9*(k-1))*sol(2+9*(k-1))*sol(4+9*(k-1))+  &
          sq2pi*sol(8+9*(k-1))*sol(5+9*(k-1))*sol(7+9*(k-1))+  &
          sq2pi*sol(8+9*(k-1))*sol(10+9*(k-1))
   bf(nline,k)=amp(k)/(std*sqrt(sol(6+9*(k-1))-sol(3+9*(k-1))+          &
    sol(4+9*(k-1))+sol(7+9*(k-1))))
!   write(0,*) "SNR: ",bf(nline,k)
enddo
do k=1,ntrace
!   dTrace(nline,k)=sol(9+9*(k-1))
   dTrace(nline,k)=sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0
enddo
!write(0,'(A2,1X,I4,3(1X,F11.3))') "**",nline,                           &
!   (sol(9+9*(k-1))+(sol(3+9*(k-1))+sol(6+9*(k-1)))/2.0d0,k=1,ntrace)
solpsf(nline,:)=sol(1:ntrace*9+1) !save PSF model

allocate(model(size(line)))
allocate(px(size(line)),py(size(line)))
do j=1,size(line)
   px(j)=real(j)  !X axis
enddo
call psfmodel1d(size(line),model,ntrace,sol)
py=real(model) !Y-axis
!   pmin=minval(py)
!   py=log10(py-pmin+1.0d0)
!cp call pgsci(6)  !change plotting colour
!cp call pgline(size(line),px,py) !plot a line
!cp call pgsci(1) !change plotting colour back to default
deallocate(px,py) !de-allocate plotting variables
deallocate(model)
!write(0,*) "pause.."
!read(5,*)

!do i=1,size(sol)
!   write(0,*) "sol:",i,sol(i)
!enddo
!read(5,*)

!fixing the shape of the PSF model (only central position and amplitude)
do k=1,ntrace
   isol(2+9*(k-1))=1 !amplitude
   isol(3+9*(k-1))=1 !position
   isol(4+9*(k-1))=1 !width
   isol(5+9*(k-1))=1 !amplitude
   isol(6+9*(k-1))=1 !position
   isol(7+9*(k-1))=1 !width
   isol(10+9*(k-1))=0 !width
enddo

allocate(solfirst(ntrace*9+1),isolfirst(ntrace*9+1))
solfirst=sol !save solution from first line
isolfirst=isol

!Next part is to use the provided positions and develop the trace.

i=nline+1
do while (i.le.naxes(1))
   line=Image(i,:)  !get next column from image to use.
   line=line-sky
   solnew=sol !save current solution
   call modelline(naxes(2),line,ntrace,solnew,isol,ilinkpsf,i) !model
!  now check the amplitudes and changes in trace
   do k=1,ntrace
      amp(k)=sq2pi*solnew(8+9*(k-1))*solnew(2+9*(k-1))*solnew(4+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(5+9*(k-1))*solnew(7+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(10+9*(k-1))
      bf(i,k)=amp(k)/max(1.0d0,(std*sqrt(solnew(6+9*(k-1))-solnew(3+9*(k-1))+ &
       solnew(4+9*(k-1))+solnew(7+9*(k-1)))))
!      if(bf(i,k).gt.0.0) write(0,*) i,k,bf(i,k)

      if((bf(i,k).lt.bcut).and.(bf(i,k).gt.0.0d0))then  !if amplitude is too low

         !try averaging together a few columns
         nbin=max(4,int((2.0d0*bcut/bf(i,k))**2.0d0)+1)
         line=0.0d0
         do j=max(i-nbin,1),min(i+nbin,naxes(1))
            line=line+Image(j,:) !Sum together ajacent lines
         enddo
         line=line/dble(2*nbin+1)-sky
         solbin=sol !try to fit again with previous solution.
         call modelline(naxes(2),line,ntrace,solbin,isol,ilinkpsf,i) !model
         amp(k)=sq2pi*solbin(8+9*(k-1))*solbin(2+9*(k-1))*solbin(4+9*(k-1))+  &
          sq2pi*solbin(8+9*(k-1))*solbin(5+9*(k-1))*solbin(7+9*(k-1))+  &
          sq2pi*solbin(8+9*(k-1))*solbin(10+9*(k-1))
         bf(i,k)=amp(k)/max(1.0d0,(std*sqrt(solbin(6+9*(k-1))-solbin(3+9*(k-1))+ &
          solbin(4+9*(k-1))+solbin(7+9*(k-1)))))

         if(bf(i,k).gt.0.0) then
            write(0,*) "binned: ",i,k,bf(i,k)
            do j=2+9*(k-1),10+9*(k-1)  !update our running solution.
               solnew(j)=solbin(j)
            enddo

         else

            do j=2+9*(k-1),10+9*(k-1)  !kill trace
               solnew(j)=0.0d0
               isol(j)=0
               bf(i,k)=0.0d0
            enddo

         endif
      endif
      tracetest=solnew(9+9*(k-1))+(solnew(3+9*(k-1))+solnew(6+9*(k-1)))/2.0d0
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
   solpsf(i,:)=sol(1:ntrace*9+1) !save PSF model for output

!   write(0,'(A2,1X,I4,3(1X,F11.3))') "**",i,(dTrace(i,k),k=1,ntrace)
   i=i+1
enddo

!restore solution from nline.
sol=solfirst
isol=isolfirst
i=nline-1
do while (i.ge.1)
   line=Image(i,:)
   line=line-sky
   solnew=sol !save current solution
   call modelline(naxes(2),line,ntrace,solnew,isol,ilinkpsf,i) !model
!  now check the amplitudes and changes in trace
   do k=1,ntrace
      amp(k)=sq2pi*solnew(8+9*(k-1))*solnew(2+9*(k-1))*solnew(4+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(5+9*(k-1))*solnew(7+9*(k-1))+  &
             sq2pi*solnew(8+9*(k-1))*solnew(10+9*(k-1))
      bf(i,k)=amp(k)/max(1.0d0,(std*sqrt(solnew(6+9*(k-1))-solnew(3+9*(k-1))+ &
       solnew(4+9*(k-1))+solnew(7+9*(k-1)))))
      if((bf(i,k).lt.bcut).and.(bf(i,k).gt.0.0d0))then  !amplitude is too low.. try binning

         !try averaging together a few columns
         nbin=max(4,int((2.0d0*bcut/bf(i,k))**2.0d0)+1)
!         write(0,*) "nbin: ",nbin
         line=0.0d0
         do j=max(i-nbin,1),min(i+nbin,naxes(1))
            line=line+Image(j,:) !Sum together ajacent lines
         enddo
         line=line/dble(2*nbin+1)-sky
         solbin=sol !try to fit again with previous solution.
         call modelline(naxes(2),line,ntrace,solbin,isol,ilinkpsf,i) !model
         amp(k)=sq2pi*solbin(8+9*(k-1))*solbin(2+9*(k-1))*solbin(4+9*(k-1))+  &
          sq2pi*solbin(8+9*(k-1))*solbin(5+9*(k-1))*solbin(7+9*(k-1))+  &
          sq2pi*solbin(8+9*(k-1))*solbin(10+9*(k-1))
         bf(i,k)=amp(k)/max(1.0d0,(std*sqrt(solbin(6+9*(k-1))-solbin(3+9*(k-1))+ &
          solbin(4+9*(k-1))+solbin(7+9*(k-1)))))

         if(bf(i,k).gt.0.0) then
            write(0,*) "binned: ",i,k,bf(i,k)
            do j=2+9*(k-1),10+9*(k-1)  !update our running solution.
               solnew(j)=solbin(j)
            enddo

         else

            do j=2+9*(k-1),10+9*(k-1)  !kill trace
               solnew(j)=0.0d0
               isol(j)=0
               bf(i,k)=0.0d0
            enddo

         endif

      endif
      tracetest=solnew(9+9*(k-1))+(solnew(3+9*(k-1))+solnew(6+9*(k-1)))/2.0d0
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
   solpsf(i,:)=sol(1:ntrace*9+1) !save PSF model for output

!   write(0,'(A2,1X,I4,3(1X,F11.3))') "**",i,(dTrace(i,k),k=1,ntrace)
   i=i-1
enddo


write(0,*) "End of Trace.. "
return
end subroutine trace
