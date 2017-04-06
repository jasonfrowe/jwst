subroutine fittransitmodel8(npars,sol,solerr)
use precision
implicit none
!import vars
integer :: npars
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr
!local vars
integer :: n,i
real(double) :: tol
real(double), allocatable, dimension(:) :: solin

allocate(solin(npars))
n=0
do i=1,npars
   if(solerr(i,1).ne.0.0)then
      n=n+1
      solin(n)=sol(i) !contains subset of sol that is passed to ldif
   endif
enddo

tol=1.0d-8 !tolerance parameter for fits

return
end subroutine fittransitmodel8
