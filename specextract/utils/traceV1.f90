subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf)
!Jason Rowe 2015 - jasonfrowe@gmail.com
!generates trace
use precision
implicit none
integer :: nkeys,nkeysmax,nKsize,i,nline,nksd2,xm,xp,j,nmaxf,k,nTrace,  &
   nTp,jj,ncut,ncutpsf,ncutd2
integer, dimension(2) :: naxes,knaxes
real(double) :: Kmin,Kmax,bpix,maxf,S,Sxy,Sxx,Sx,Sy,df,bcut,f
real(double), dimension(:,:) :: Image,dTrace,bf
real(double), allocatable, dimension(:) :: line,lpsf,lpsftemp,psfwork
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

!constants
nKsize=64        !size of Kernel for import (note.. this gets resized)
nksd2=nKsize/2 !precompute size of PSF/2
!parameters to control trace
ncut=35          !width of spectrum to zero out
bcut=1.0d0    !threshold for finding a trace
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

!allocate(lpsftemp(size(lpsf(nKsize/2-ncutpsf:nKsize/2+ncutpsf))))
!lpsftemp=lpsf(nKsize/2-ncutpsf:nKsize/2+ncutpsf)
!nKsize=size(lpsf(nKsize/2-ncutpsf:nKsize/2+ncutpsf))
!deallocate(lpsf)
!call move_alloc(lpsftemp,lpsf) !move temp to psf and deallocate temp
!write(0,*) "Size: ",size(lpsf)

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
!   if(k.eq.1)then
!      do j=1,size(line(xm:xp))
!         write(6,*) px(j),py(j)
!      enddo
!   endif
   deallocate(px,py) !de-allocate plotting variables


!There can be siginficant differences compare to first pass, so
!zero out residuals to avoid problems.
   xm=max(1,i-ncut+1)
   xp=min(naxes(2),i+ncut)
   do j=xm,xp
      line(j)=0.0d0 !zero out value - should actually be set to sky.
   enddo
enddo


!Update the PSF. Start by getting the orignal line from the Image
line=Image(nline,:)    !we start at user-defined 'nline'
!subtract off minimum value - maybe a true sky value would be better
line=line-minval(line)
psfwork=0.0d0 !initialize work
do k=1,ntrace !look over k
   i=dTrace(nline,k)
   xm=max(1,i-nksd2+1)      !we now use the data to update the model PSF
   xp=min(naxes(2),i+nksd2)
   psfwork(nksd2-(i-xm):nksd2+(xp-i))=psfwork(nksd2-(i-xm):nksd2+(xp-i))+ &
      line(xm:xp)*bf(nline,k)
enddo
psfwork=psfwork/sum(bf(nline,1:ntrace))
do i=1,nTrace    !make a trace for each order.
   psf(:,i)=psfwork(:)  !our initial guess for each order, now refined
enddo

!Next part is to use the found positions and develop the trace.

do i=nline+1,naxes(1) !forward direction
   line=Image(i,:)  !get next column from image to use.
   line=line-minval(line) !offset to zero - change to sky-subtraction

   do k=1,nTrace
      nTp=dTrace(i-1,k)
!      write(0,*) "*",k,nTp
      maxf=0.0d0
      nmaxf=dTrace(i-1,k)
      dTrace(i,k)=dTrace(i-1,k) !pre-set
      do j=max(1,nTp-1),min(nTp+1,naxes(2))
         xm=max(1,j-nksd2+1)
         xp=min(naxes(2),j+nksd2)
         df=Sum(line(xm:xp)*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
         if(df.gt.maxf)then
            maxf=df
            nmaxf=j
         endif
      enddo
      nTp=nmaxf
!      write(0,*) "next."
      !position of trace
      j=nmaxf
      xm=max(1,j-nksd2+1)
      xp=min(naxes(2),j+nksd2)
      !calculate best-fit scaling
      S=1.0d0
      Sxy=Sum(line(xm:xp)*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sxx=Sum(lpsf(nksd2-(j-xm):nksd2+(xp-j))*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sx =Sum(lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sy =Sum(line(xm:xp))
      bf(i,k)=(S*Sxy-Sx*Sy)/(S*Sxx-Sx*Sx)
!      write(0,*) "** next",bf
      if(bf(i,k).gt.bcut)then
         dTrace(i,k)=nTp !save next entry from trace
!      write(0,*) "bf: ",k,nTp,bf
!         do jj=xm,xp
!            line(jj)=line(jj)-lpsf(nksd2-(i-jj))*bf
!         enddo
         line(xm:xp)=line(xm:xp)-lpsf(nksd2-(j-xm):nksd2+(xp-j))*bf(i,k)
!zero out residuals
         xm=max(1,j-ncut+1)
         xp=min(naxes(2),j+ncut)
         do jj=xm,xp
            line(jj)=0.0!line(j)-lpsf(nksd2-(i-j))*bf
         enddo
      else
         dTrace(i,k)=dTrace(i-1,k)!0.0d0
      endif

   enddo
!   write(0,*) "**",i,(dTrace(i,k),k=1,3)

enddo

!nTp=nTstart(k)
do i=nline-1,1,-1 !negative direction
   line=Image(i,:)
   line=line-minval(line)

   do k=1,nTrace
      nTp=dTrace(i+1,k)
!      write(0,*) "k:",k,nTp
      maxf=0.0d0
      nmaxf=dTrace(i+1,k)
      dTrace(i,k)=dTrace(i+1,k) !pre-set
      do j=max(1,nTp-1),min(nTp+1,naxes(2))
         xm=max(1,j-nksd2+1)
         xp=min(naxes(2),j+nksd2)
         df=Sum(line(xm:xp)*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
         if(df.gt.maxf)then
            maxf=df
            nmaxf=j
         endif
      enddo
      nTp=nmaxf

      !position of trace
      j=nmaxf
      xm=max(1,j-nksd2+1)
      xp=min(naxes(2),j+nksd2)
      !calculate best-fit scaling
      S=1.0d0
      Sxy=Sum(line(xm:xp)*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sxx=Sum(lpsf(nksd2-(j-xm):nksd2+(xp-j))*lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sx =Sum(lpsf(nksd2-(j-xm):nksd2+(xp-j)))
      Sy =Sum(line(xm:xp))
      bf(i,k)=(S*Sxy-Sx*Sy)/(S*Sxx-Sx*Sx)

      if(bf(i,k).gt.bcut)then
         dTrace(i,k)=nTp !save next entry from trace
!      write(0,*) "bf: ",k,nTp,bf
!         do jj=xm,xp
!            line(jj)=line(jj)-lpsf(nksd2-(i-jj))*bf
!         enddo
         line(xm:xp)=line(xm:xp)-lpsf(nksd2-(j-xm):nksd2+(xp-j))*bf(i,k)
!zero out residuals
         xm=max(1,j-ncut+1)
         xp=min(naxes(2),j+ncut)
         do jj=xm,xp
            line(jj)=0.0!line(j)-lpsf(nksd2-(i-j))*bf
         enddo
      else
         dTrace(i,k)=dTrace(i+1,k)!0.0d0
      endif

   enddo
!   write(0,*) i,(dTrace(i,k),k=1,3)

enddo

!do i=1,naxes(1)
!   write(6,*) i,dTrace(i,k)
!enddo

!plotting to have a look
!allocate(px(naxes(2)),py(naxes(2)))
!do i=1,naxes(2)
!   px(i)=real(i)
!enddo
!py=real(f)
!!pmin=minval(py)
!!py=log10(py-pmin+1.000)
!call pgpage()
!call pgsci(1)
!call pgvport(0.10,0.95,0.15,0.95) !make room around the edges for labels
!call pgwindow(minval(px),maxval(px),minval(py),maxval(py)) !plot scale
!call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
!call pglabel("X (pixels)","Y (Counts)","")
!call pgline(naxes(2),px,py)
!deallocate(px,py)

write(0,*) "End of Trace.. "
return
end subroutine trace
