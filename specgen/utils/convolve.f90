!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine convolve(xmax,ymax,pixels,nrK,nKs,rKernel,noversample,       &
   cpixels,ybounds,ntrace)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: xmax,ymax,nK,i,j,l,m,nKd2,ii,jj,ii1,ii2,jj1,jj2,nrK,nKs,     &
   noversample,Ks,ntrace
integer, dimension(2) :: ybounds
real(double) :: mSum,p2w,wl,p
real(double), dimension(xmax,ymax) :: pixels,cpixels
real(double), allocatable, dimension(:,:) :: subIm,subIm2,Kernel,wKernel
real(double), dimension(nrK,nKs,nKs) :: rKernel
character(80) :: line

interface
   subroutine genkernel(nrK,nKs,rKernel,Kernel,wl)
      use precision
      integer, intent(inout) :: nrK,nKs
      real(double), intent(inout) :: wl
      real(double), dimension(:,:),intent(inout) :: Kernel
      real(double), dimension(:,:,:), intent(inout) :: rKernel
   end subroutine
end interface

!write(0,*) "Convolution begins"

Ks=50*noversample !this stub allows for a smaller piece of the Kernel to be used

!adjust ybounds for Kernel width for convolution speed up
nKd2=Ks/2 !half size of kernel
ybounds(1)=max(1,ybounds(1)-nKd2)
ybounds(2)=min(ymax,ybounds(2)+nKd2)
write(0,*) 'ybounds: ',ybounds(1),ybounds(2)

cpixels=0.0d0 !initialize convoluted image to zero
allocate(subIm(Ks,Ks),subIm2(Ks,Ks)) !sub-Image stamps for multiplication

do i=1,xmax

!  This part is where we generate the wavelength specific Kernel
   if(mod(i,1).eq.0)then !only need to update Kernel every 10th pixel
      nK=nKs !native Kernel should already match image scale
      if(allocated(Kernel)) deallocate(Kernel)
      allocate(Kernel(nK,nK))
      Kernel=0.0d0
      p=dble(i)
      wl=p2w(p,noversample,ntrace)/10000.0d0 !A -> um
!      write(0,*) "new Kernel wl: ",wl
      call genkernel(nrK,nKs,rKernel,Kernel,wl)
      allocate(wKernel(Ks,Ks))
      wKernel=Kernel(nK/2-Ks/2:nK/2-Ks/2+Ks,nK/2-Ks/2:nK/2-Ks/2+Ks)
!      write(0,*) "Ksize: ",size(wKernel,1),size(wKernel,2)
      deallocate(Kernel)
      allocate(Kernel(Ks,Ks))
      nK=Ks
      Kernel=wKernel
      deallocate(wKernel)
      mSum=Sum(Kernel)
      Kernel=Kernel/mSum
   endif

   write(line,502) "Convolution: ",xmax,i,wl
   502 format(A13,I5,1X,I5,1X,F6.3)
   call ovrwrt(line,2)
!   write(0,'(A80)') line
   do j=ybounds(1),ybounds(2) !only scan area with non-zero pixels
      subIm=0.0d0 !initialize subImage
      ii=0 !ii,jj mark sub-Image position
      do l=i-nkd2,i-nkd2+nK-1
         ii=ii+1
         jj=0
         do m=j-nkd2,j-nkd2+nK-1
            jj=jj+1
            !check that pixel location is valid
            if((l.gt.0).and.(l.le.xmax).and.(m.gt.0).and.(m.le.ymax))then
!               write(0,*) ii,jj,l,m,nK
               subIm(ii,jj)=pixels(l,m) !copy image value to stamp
            endif
         enddo
      enddo
      mSum=Sum(abs(subIm))

      if(mSum.gt.0.0e-30)then !if array is full of zeros, skip..
!      write(0,*) "multiply ",i,j
         subIm2=matmul(subIm,Kernel) !multiple subImage by Kernel
         ii=0
         do l=i-nkd2,i-nkd2+nK-1
            ii=ii+1
            jj=0
            do m=j-nkd2,j-nkd2+nK-1
               jj=jj+1
               if((l.gt.0).and.(l.le.xmax).and.(m.gt.0).and.(m.le.ymax))then
                  cpixels(l,m)=cpixels(l,m)+subIm2(ii,jj) !copy convolved stamp
               endif
            enddo
         enddo
      endif
   enddo
enddo
write(0,*) " " !clear next line of text output from ovrwrt usage above

return
end
