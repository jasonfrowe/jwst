subroutine trace(naxes,Image,bpix,nline,nTrace,dTrace,bf)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nkeys,nkeysmax,nKsize,i,nline,nksd2,xm,xp,j,nmaxf,k,nTrace,  &
   nTp,jj,ncut
integer, allocatable, dimension(:) :: nTstart
integer, dimension(2) :: naxes,knaxes
real(double) :: Kmin,Kmax,bpix,maxf,S,Sxy,Sxx,Sx,Sy,df,bcut
real(double), dimension(:,:) :: Image,dTrace,bf
real(double), allocatable, dimension(:) :: lpsf,line,f
real(double), allocatable, dimension(:,:) :: Kernel
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

ncut=25           !width of spectrum to zero out
bcut=10000.0d0    !threshold for finding a trace

!read in a Kernel
kfile="Kernels1/psf_2.10000_m1_dx0.69_dy0.32_LLNLCoated.fits"
nkeysmax=800
allocate(header(nkeysmax))
nKsize=256
allocate(Kernel(nKsize,nKsize))
call getfits(kfile,knaxes,Kernel,Kmin,Kmax,nkeys,header,bpix)
write(0,*) "knaxes: ",knaxes

!generate 1-D PSF for cross-correlation
allocate(lpsf(nKsize))
lpsf=Sum(Kernel,2)
lpsf=lpsf-minval(lpsf)

!allocate space for line-scan
allocate(line(naxes(2)))
nksd2=nKsize/2 !size of PSF/2

line=Image(nline,:)
line=line-minval(line)
allocate(f(naxes(2)))

allocate(nTstart(nTrace))
do k=1,nTrace !loop over expected number of traces

   maxf=0.0d0
   do i=1,naxes(2) !scan over the line
      xm=max(1,i-nksd2+1)
      xp=min(naxes(2),i+nksd2)
      f(i)=Sum(line(xm:xp)*lpsf(nksd2-(i-xm):nksd2+(xp-i)))

      if(f(i).gt.maxf)then
         maxf=f(i)
         nmaxf=i
      endif
   enddo
   write(0,*) "maxf: ",nmaxf,maxf
   nTstart(k)=nmaxf !initial position of trace

   !position of trace
   i=nmaxf
   xm=max(1,i-nksd2+1)
   xp=min(naxes(2),i+nksd2)
   !calculate best-fit scaling
   S=1.0d0
   Sxy=Sum(line(xm:xp)*lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sxx=Sum(lpsf(nksd2-(i-xm):nksd2+(xp-i))*lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sx =Sum(lpsf(nksd2-(i-xm):nksd2+(xp-i)))
   Sy =Sum(line(xm:xp))
   bf(nline,k)=(S*Sxy-Sx*Sy)/(S*Sxx-Sx*Sx)
!   write(0,*) "bf: ",bf
!   do j=xm,xp
!      line(j)=line(j)-lpsf(nksd2-(i-j))*bf
!   enddo
   line(xm:xp)=line(xm:xp)-lpsf(nksd2-(i-xm):nksd2+(xp-i))*bf(nline,k)
!   write(0,*) "Pass.."

!zero out residuals
   xm=max(1,i-ncut+1)
   xp=min(naxes(2),i+ncut)
   do j=xm,xp
      line(j)=0.0!line(j)-lpsf(nksd2-(i-j))*bf
   enddo

enddo
deallocate(f)

do k=1,nTrace
   dTrace(nline,k)=nTstart(k)
enddo
!write(0,*) "nTstart: ",nTstart(k)
!nTp=nTstart(k)
do i=nline+1,naxes(1) !forward direction
   line=Image(i,:)
   line=line-minval(line) !offset to zero.

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
