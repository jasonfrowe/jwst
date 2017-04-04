!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getnumfitpars(nunit,npars)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!input vars
integer :: nunit,npars
!local vars
integer :: i,nvar,filestatus,iltype
real(double) :: dumr
character(3) :: pname

npars=0 !counting number of model parameters
i=0 !counter for number of observations.
iltype=0 !if =0, then line describes variable, otherwise it's parameter values
do
   if (iltype.eq.0) then
      read(nunit,*,iostat=filestatus) pname,nvar
   else
      read(nunit,*,iostat=filestatus) dumr
   endif
   if(filestatus == 0) then
      i=i+1
      if(iltype.eq.0)then
         npars=npars+nvar
         iltype=nvar
      else
         iltype=iltype-1
      endif
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

rewind(nunit)

return
end subroutine getnumfitpars

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getfitpars()
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none

write(0,*) "hello"

return
end subroutine getfitpars
