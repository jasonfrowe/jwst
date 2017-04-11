subroutine exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
use precision
implicit none
!import vars
integer :: nunit,npars,nplanet
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr
!local vars
integer :: i,j,k,np,ii

do i=1,8
   select case (i)
      case (1)
         write(nunit,500) "RHO",solrange(i,2)-solrange(i,1)+1
      case (2:5)
         write(nunit,501) "NL",i-1,solrange(i,2)-solrange(i,1)+1
      case (6)
         write(nunit,500) "DIL",solrange(i,2)-solrange(i,1)+1
      case (7)
         write(nunit,500) "VOF",solrange(i,2)-solrange(i,1)+1
      case (8)
         write(nunit,500) "ZPT",solrange(i,2)-solrange(i,1)+1
   end select
   do j=solrange(i,1),solrange(i,2)
      write(nunit,502) sol(j),(solerr(j,k),k=1,4)
   enddo
enddo
500 format(A3,1X,I6)
501 format(A2,I1,1X,I6)
502 format(5(1PE17.10,1X))

!planet parameters
do np=1,nplanet !loop over each planet in transit model
   do i=1,10    !loop over each parameter for a single planet
      ii=8+i+nplanet*(np-1)
      select case (i)
         case (1)
            write(nunit,501) "EP",np,solrange(ii,2)-solrange(ii,1)+1
         case (2)
            write(nunit,501) "PE",np,solrange(ii,2)-solrange(ii,1)+1
         case (3)
            write(nunit,501) "BB",np,solrange(ii,2)-solrange(ii,1)+1
         case (4)
            write(nunit,501) "RD",np,solrange(ii,2)-solrange(ii,1)+1
         case (5)
            write(nunit,501) "EC",np,solrange(ii,2)-solrange(ii,1)+1
         case (6)
            write(nunit,501) "ES",np,solrange(ii,2)-solrange(ii,1)+1
         case (7)
            write(nunit,501) "KR",np,solrange(ii,2)-solrange(ii,1)+1
         case (8)
            write(nunit,501) "TE",np,solrange(ii,2)-solrange(ii,1)+1
         case (9)
            write(nunit,501) "EL",np,solrange(ii,2)-solrange(ii,1)+1
         case (10)
            write(nunit,501) "AL",np,solrange(ii,2)-solrange(ii,1)+1
      end select
      do j=solrange(ii,1),solrange(ii,2)
         write(nunit,502) sol(j),(solerr(j,k),k=1,4)
      enddo
   enddo
enddo

return
end subroutine exportfitpars
