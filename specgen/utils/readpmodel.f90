subroutine readpmodel(nunit,nmodel,rprs)
use precision
implicit none
integer :: nunit,nmodel
real(double), dimension(:) :: rprs
!local vars
integer :: npt,i,filestatus
real(double), allocatable, dimension(:) :: wv,pmod

rprs=0.1188

allocate(wv(nmodel),pmod(nmodel))
!i=1
do
   if(i.gt.nmodel)then !ran out of array space
      write(0,*) "Critical Error: Planet model has higher resolution that Star model"
      stop
   endif
   read(nunit,*) wv(i),pmod(i)
   if(filestatus == 0) then
      i=i+1
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

return
end subroutine readpmodel
