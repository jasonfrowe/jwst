program genoskernel
use precision
implicit none
integer noversample,nrK,nKs,i,j,k,ii,jj,l,m,nK
real(double) :: dnover,x1,x2,y
real(double), allocatable, dimension(:) :: x1a,x2a
real(double), allocatable, dimension(:,:) :: wKernel,Kernel,y2a
real(double), allocatable, dimension(:,:,:) :: rKernel
character(80), allocatable, dimension(:) :: filenames

interface
   subroutine readKernels(nrK,nK,rKernel,noversample)
      use precision
      implicit none
      integer,intent(in) :: noversample
      integer,intent(inout) :: nrK,nK
      real(double), dimension(:,:,:), intent(inout) :: rKernel
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

!oversampling
noversample=8
dnover=dble(noversample) !double int -> double

!read in Kernels
nrK=26 !number of Kernels to readin
nKs=256 !native size of Kernels
allocate(rKernel(nrK,nKs,nKs))
call readKernels(nrK,nKs,rKernel,1) !read in native size Kernels

!get filenames for output
allocate(filenames(nrK))
call getfilenames(noversample,nrK,filenames)

!allocate space for oversample Kernels
allocate(wKernel(nKs,nKs))
nK=nKs*noversample
allocate(Kernel(nK,nK))

!set up arrays for spline calculations
allocate(x1a(nKs),x2a(nKs)) !allocate array for spline X,Y co-ordinates
do i=1,nKs
   x1a(i)=dble(i) !co-ordinates for grid to spline
enddo
x2a=x1a !copy x1a to x2a, since they are the same (square-array)
allocate(y2a(nKs,nKs)) !array to hold derivatives from splie2

do k=1,nrK
   write(0,501) "Kernel #",k,"/",nrK
   write(0,502) filenames(k)
   501 format(A8,I2,A1,I2)
   502 format(A80)
   wKernel=rKernel(k,:,:)
!  use a cubic spline to oversample the Kernel
   call splie2(x1a,x2a,wKernel,nKs,nKs,y2a)  !calculate derivatives
   do i=1,nKs
      do j=1,nKs
         ii=i*noversample-noversample+1
         jj=j*noversample-noversample+1
         do l=1,noversample
            do m=1,noversample
               x1=dble(ii+l-1)/dnover
               x2=dble(jj+m-1)/dnover
               call splin2(x1a,x2a,wKernel,y2a,nKs,nKs,x1,x2,y)
               Kernel(ii+l-1,jj+m-1)=y
!            write(0,*) i,j,wKernel(i,j),Kernel(ii+l-1,jj+m-1)
!            read(5,*)
            enddo
         enddo
      enddo
   enddo
!  write the new oversampled Kernel
   call writefits(nK,nK,Kernel,filenames(k))
enddo

end program genoskernel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getfilenames(noversample,nrK,filenames)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer noversample,nrK
character(80), dimension(nrK) :: filenames
!filenames of Kernels

write(filenames(1),500) "Kernels",noversample,"/psf_0.50000_m1_dx0.52_dy0.16_LLNLCoated.fits"
write(filenames(2),500) "Kernels",noversample,"/psf_0.60000_m1_dx0.11_dy0.08_LLNLCoated.fits"
write(filenames(3),500) "Kernels",noversample,"/psf_0.70000_m1_dx0.93_dy0.88_LLNLCoated.fits"
write(filenames(4),500) "Kernels",noversample,"/psf_0.80000_m1_dx0.21_dy0.61_LLNLCoated.fits"
write(filenames(5),500) "Kernels",noversample,"/psf_0.90000_m1_dx0.10_dy0.84_LLNLCoated.fits"
write(filenames(6),500) "Kernels",noversample,"/psf_1.00000_m1_dx0.24_dy0.27_LLNLCoated.fits"
write(filenames(7),500) "Kernels",noversample,"/psf_1.10000_m1_dx0.87_dy0.53_LLNLCoated.fits"
write(filenames(8),500) "Kernels",noversample,"/psf_1.20000_m1_dx0.12_dy0.79_LLNLCoated.fits"
write(filenames(9),500) "Kernels",noversample,"/psf_1.30000_m1_dx0.49_dy0.63_LLNLCoated.fits"
write(filenames(10),500) "Kernels",noversample,"/psf_1.40000_m1_dx0.78_dy0.53_LLNLCoated.fits"
write(filenames(11),500) "Kernels",noversample,"/psf_1.50000_m1_dx0.44_dy0.38_LLNLCoated.fits"
write(filenames(12),500) "Kernels",noversample,"/psf_1.60000_m1_dx0.48_dy0.99_LLNLCoated.fits"
write(filenames(13),500) "Kernels",noversample,"/psf_1.70000_m1_dx0.17_dy0.56_LLNLCoated.fits"
write(filenames(14),500) "Kernels",noversample,"/psf_1.80000_m1_dx0.77_dy0.23_LLNLCoated.fits"
write(filenames(15),500) "Kernels",noversample,"/psf_1.90000_m1_dx0.43_dy0.08_LLNLCoated.fits"
write(filenames(16),500) "Kernels",noversample,"/psf_2.00000_m1_dx0.00_dy0.00_LLNLCoated.fits"
write(filenames(17),500) "Kernels",noversample,"/psf_2.10000_m1_dx0.69_dy0.32_LLNLCoated.fits"
write(filenames(18),500) "Kernels",noversample,"/psf_2.20000_m1_dx0.69_dy0.38_LLNLCoated.fits"
write(filenames(19),500) "Kernels",noversample,"/psf_2.30000_m1_dx0.49_dy0.58_LLNLCoated.fits"
write(filenames(20),500) "Kernels",noversample,"/psf_2.40000_m1_dx0.16_dy0.97_LLNLCoated.fits"
write(filenames(21),500) "Kernels",noversample,"/psf_2.50000_m1_dx0.77_dy0.26_LLNLCoated.fits"
write(filenames(22),500) "Kernels",noversample,"/psf_2.60000_m1_dx0.30_dy0.38_LLNLCoated.fits"
write(filenames(23),500) "Kernels",noversample,"/psf_2.70000_m1_dx0.76_dy0.36_LLNLCoated.fits"
write(filenames(24),500) "Kernels",noversample,"/psf_2.80000_m1_dx0.11_dy0.44_LLNLCoated.fits"
write(filenames(25),500) "Kernels",noversample,"/psf_2.90000_m1_dx0.36_dy0.52_LLNLCoated.fits"
write(filenames(26),500) "Kernels",noversample,"/psf_3.00000_m1_dx0.60_dy0.50_LLNLCoated.fits"
500 format(A7,I1,A45)

return
end
