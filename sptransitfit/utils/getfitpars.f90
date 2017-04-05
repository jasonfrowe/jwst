!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getnumfitpars(nunit,npars,nplanetmax)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!(c) Jason Rowe 2017
use precision
implicit none
!input vars
integer :: nunit,npars,nplanetmax
!local vars
integer :: i,nvar,filestatus,iltype,np
real(double) :: dumr
character(3) :: pname

nplanetmax=0 !counting number of planets in the model
npars=0 !counting number of model parameters
i=0 !counter for number of observations.
iltype=0 !if =0, then line describes variable, otherwise it's parameter values
do
   if (iltype.eq.0) then
      read(nunit,*,iostat=filestatus) pname,nvar
      if((pname(1:2).eq.'EP').or.(pname(1:2).eq.'PE').or.               &
       (pname(1:2).eq.'BB').or.(pname(1:2).eq.'RD').or.                 &
       (pname(1:2).eq.'EC').or.(pname(1:2).eq.'ES').or.                 &
       (pname(1:2).eq.'KR').or.(pname(1:2).eq.'TE').or.                 &
       (pname(1:2).eq.'EL').or.(pname(1:2).eq.'AL'))then
         read(pname(3:3),*) np !read in planet number
         nplanetmax=max(nplanetmax,np) !get numbe of planets in model
      endif
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
subroutine getfitpars(nunit,nparsmax,nplanetmax,npars,nplanet,sol,      &
 solerr,solrange)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!(c) Jason Rowe 2017
use precision
implicit none
!import vars
integer :: nunit,nparsmax,npars,nplanet,nplanetmax
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,solrange
!local vars
integer :: i,j,k,iltype,nvar,filestatus,np
real(double), allocatable, dimension(:) :: parsin
character(3) :: pname

allocate(parsin(5))

nplanet=0 !counting number of planets in model
npars=0 !counting number of model parameters
i=0 !counter for number of observations.
j=0 !counter for number of model parameters
iltype=0 !if =0, then line describes variable, otherwise it's parameter values
do
   if (iltype.eq.0) then
      read(nunit,*,iostat=filestatus) pname,nvar
      if((pname(1:2).eq.'EP').or.(pname(1:2).eq.'PE').or.               &
       (pname(1:2).eq.'BB').or.(pname(1:2).eq.'RD').or.                 &
       (pname(1:2).eq.'EC').or.(pname(1:2).eq.'ES').or.                 &
       (pname(1:2).eq.'KR').or.(pname(1:2).eq.'TE').or.                 &
       (pname(1:2).eq.'EL').or.(pname(1:2).eq.'AL'))then
         read(pname(3:3),*) np !read in planet number
         nplanet=max(nplanet,np) !get numbe of planets in model
      endif
   else
      read(nunit,*,iostat=filestatus) (parsin(k),k=1,5)
   endif
   if(filestatus == 0) then
      i=i+1
      if(iltype.eq.0)then
         npars=npars+nvar
         iltype=nvar
      else
         j=j+1
         if(j.gt.nparsmax)then !make sure we are not overunning array
            write(0,*) "nparsmax is too small. This should not happen!"
            write(0,*) "nparsmax: ",nparsmax,j
            stop
         else
            sol(j)=parsin(1)
            solerr(j,:)=parsin(2:5)
            iltype=iltype-1
         endif
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

return
end subroutine getfitpars
