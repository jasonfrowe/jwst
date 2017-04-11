subroutine exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
use precision
implicit none
!import vars
integer :: nunit,npars,nplanet
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr
!local vars
integer :: i

do i=1,8
   select case (i)
      case (1)
         write(nunit,500) "RHO",solrange(i,2)-solrange(i,1)+1
      case (2:5)
         write(nunit,501) "NL",i-1,solrange(i,2)-solrange(i,1)+1
   end select
enddo
500 format(A3,1X,I6)
501 format(A2,I1,1X,I6)

return
end subroutine exportfitpars
