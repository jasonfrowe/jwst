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
   else
      read(nunit,*,iostat=filestatus) dumr
   endif
   if(filestatus == 0) then
      i=i+1
      if(iltype.eq.0)then
         if((pname(1:2).eq.'EP').or.(pname(1:2).eq.'PE').or.            &
          (pname(1:2).eq.'BB').or.(pname(1:2).eq.'RD').or.              &
          (pname(1:2).eq.'EC').or.(pname(1:2).eq.'ES').or.              &
          (pname(1:2).eq.'KR').or.(pname(1:2).eq.'TE').or.              &
          (pname(1:2).eq.'EL').or.(pname(1:2).eq.'AL'))then
            read(pname(3:3),*) np !read in planet number
            nplanetmax=max(nplanetmax,np) !get numbe of planets in model
         endif
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
integer :: i,j,k,iltype,nvar,filestatus,np,ipar,getparnum
real(double), allocatable, dimension(:) :: parsin
character(3) :: pname

allocate(parsin(5))

solrange=0 !init to zero.  If indices are zero, then model parameter is
           !undefined
nplanet=0 !counting number of planets in model
npars=0 !counting number of model parameters
i=0 !counter for number of observations.
j=0 !counter for number of model parameters
ipar=0 !counter for base of model parameters (RHO,NL1,NL2,...)
iltype=0 !if =0, then line describes variable, otherwise it's parameter values
do
   if (iltype.eq.0) then
      read(nunit,*,iostat=filestatus) pname,nvar
   else
      read(nunit,*,iostat=filestatus) (parsin(k),k=1,5)
   endif
   if(filestatus == 0) then
      i=i+1
      if(iltype.eq.0)then
         !the getparnum function compares pname to list of parameters and
         !returns appropriate index for that model parameter
         !source is found in this file
         ipar=getparnum(pname)
         write(0,'(A3,1X,I3,1X,I5,1X,I5,1X,I5)') pname,ipar,nvar,       &
          npars+1,npars+nvar
         solrange(ipar,1)=npars+1
         solrange(ipar,2)=npars+nvar
         if((pname(1:2).eq.'EP').or.(pname(1:2).eq.'PE').or.            &
          (pname(1:2).eq.'BB').or.(pname(1:2).eq.'RD').or.              &
          (pname(1:2).eq.'EC').or.(pname(1:2).eq.'ES').or.              &
          (pname(1:2).eq.'KR').or.(pname(1:2).eq.'TE').or.              &
          (pname(1:2).eq.'EL').or.(pname(1:2).eq.'AL'))then
            read(pname(3:3),*) np !read in planet number
            nplanet=max(nplanet,np) !get numbe of planets in model
         endif
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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
integer function getparnum(pname)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!scan though parameter names and get proper indice.
implicit none
!import vars
character(3) :: pname
!local vars
integer :: i,np
character(3), dimension(18) :: names
data names /'RHO','NL1','NL2','NL3','NL4','DIL','VOF','ZPT','EPO','PER',&
'BBB','RDR','ECW','ESW','KRV','TED','ELL','ALB'/

getparnum=0 !default to zero.  If zero is returned then we have a problem

do i=1,8
   if(names(i).eq.pname)then
      getparnum=i
      return
   endif
enddo

do i=9,18
   if(names(i)(1:2).eq.pname(1:2))then
      read(pname(3:3),*) np
      getparnum=i+(np-1)*10
      return
   endif
enddo

return
end
