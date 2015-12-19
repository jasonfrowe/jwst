function getsky(naxes,Image,std)
use precision
implicit none
integer i,j,nmax,status,k,npt,itmax,ncut
integer, dimension(2) :: naxes
integer, allocatable, dimension(:) :: ip
real(double) :: getsky,median,std,stdev2,stdcut
real(double), allocatable, dimension(:) :: pixels
real(double), dimension(:,:) :: Image

stdcut=3.0

getsky=0.0d0

nmax=naxes(1)*naxes(2)
allocate(pixels(nmax),stat=status)
if(status.gt.0) then !check that large array was allocated
   write(0,*) "Allocation of Image array failed.."
   write(0,*) "Status: ",status
   stop
endif

!create vector with the pixel values
k=0
do i=1,naxes(1)
   do j=1,naxes(2)
      k=k+1
      pixels(k)=Image(i,j)
   enddo
enddo
npt=k

!sort data to get median
allocate(ip(npt))
call rqsort(npt,pixels,ip)
median=pixels(ip(npt/2))
!write(0,*) "Median: ",median

!calcute standard-deviation around the median
std=stdev2(npt,pixels,median)
!write(0,*) "std: ",std


itmax=30
j=1
ncut=1
do while((ncut.gt.0).and.(j.lt.itmax))

   !sigma-clip data
   k=0
   do i=1,npt
      if(abs(pixels(i)-median).lt.stdcut*std)then
         k=k+1
         pixels(k)=pixels(i)
      endif
   enddo
   ncut=npt-k
   npt=k

   !recalculate median
   call rqsort(npt,pixels,ip)
   median=pixels(ip(npt/2))
!   write(0,*) "Median: ",median

   !calcute standard-deviation around the median
   std=stdev2(npt,pixels,median)
!   write(0,*) "std: ",std,j,ncut

!   write(0,*) "Skycal ",median,std
   j=j+1 !increase counter
enddo
getsky=median

!read(5,*)

return
end function getsky
