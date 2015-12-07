subroutine readKernels(nrK,nK,rKernel,noversample)
use precision
implicit none
integer :: nrK,nK,nkeys,nkeysmax,i,j,noversample
integer, dimension(2) :: naxes
real(double) :: Krmin,Krmax,bpix
real(double), allocatable, dimension(:,:) :: Kernel
real(double), dimension(:,:,:) :: rKernel
character(80), allocatable, dimension(:) :: filenames,header

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

bpix=1000000.0 !marking bad pixels

!filenames of Kernels
allocate(filenames(nrK))
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
!there seems to be an error with the 2 um PSF.. so replacing it with 2.1 for now
!write(filenames(16),500) "Kernels",noversample,"/psf_2.00000_m1_dx0.00_dy0.00_LLNLCoated.fits"
write(filenames(16),500) "Kernels",noversample,"/psf_2.10000_m1_dx0.69_dy0.32_LLNLCoated.fits"
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

nkeysmax=700
allocate(Kernel(nK,nK),header(nkeysmax))

do i=1,nrK
   call getfits(filenames(i),naxes,Kernel,Krmin,Krmax,nkeys,header,bpix)
!   write(0,*) i,Krmin,Krmax
   rKernel(i,:,:)=transpose(Kernel(:,:))
enddo


return
end
