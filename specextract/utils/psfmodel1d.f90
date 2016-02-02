subroutine psfmodel1d(npt,model,ntrace,sol)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: npt,ntrace,k,nfit,i
real(double) x,triplegaussian
real(double), dimension(:) :: model,sol
real(double), allocatable, dimension(:) :: sol1

!write(0,*) "PSFModel?? start"

nfit=10 !parameters for triple-gaussian
allocate(sol1(nfit))

model=0 !initalize output array for model values to zero
sol1(1)=sol(1) !SKY
do k=1,ntrace
   sol1(2:10)=sol(2+9*(k-1):10+9*(k-1))
   do i=1,npt
      x=dble(i)
      model(i)=model(i)+triplegaussian(nfit,sol1,x)
!      write(6,*) x,model(i)
   enddo
enddo

!write(0,*) "PSFModel?? end"

end subroutine psfmodel1d
