subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,     &
 exptime,sptmodel)
use precision
implicit none
!import vars
integer :: nplanet,npars,nwv,nobs
real(double), dimension(:) :: sol,time,exptime
real(double), dimension(:,:) :: solrange,sptmodel
!local vars
integer :: nintg,iwv
real(double) :: Pi,tPi,pid2,G,Cs,fDB,c1,c2,c3,c4,dil,voff,zpt

interface
   subroutine getbasicpars(iwv,sol,solrange,c1,c2,c3,c4,dil,voff,zpt)
      use precision
      implicit none
      integer :: iwv
      real(double) :: c1,c2,c3,c4,dil,voff,zpt
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solrange
   end subroutine getbasicpars
end interface

!Model parameters
nintg=41 !number of samples to convolve integration time

!Physical Constants
Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi
pid2=Pi/2.0d0  !and Pi/2
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Cs=2.99792458e8 !Speed of light
fDB=1.0 !Doppler Boosting factor

do iwv=1,nwv !loop over all bandpasses

   call getbasicpars(iwv,sol,solrange,c1,c2,c3,c4,dil,voff,zpt)



   read(5,*) !simple pause statement

enddo


return
end subroutine sptransitmodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getbasicpars(iwv,sol,solrange,c1,c2,c3,c4,dil,voff,zpt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: iwv
real(double) :: c1,c2,c3,c4,dil,voff,zpt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solrange
!local vars
integer i
integer, parameter :: npars=7
real(double), dimension(npars) :: lpars

do i=1,npars
   if(solrange(i+1,2)-solrange(i+1,1).gt.0)then
      lpars(i)=sol(solrange(i+1,1)+iwv-1)
   else
      lpars(i)=sol(solrange(i+1,1))
   endif
enddo

c1=lpars(1)
c2=lpars(2)
c3=lpars(3)
c4=lpars(4)
dil=lpars(5)
voff=lpars(6)
zpt=lpars(7)

return
end
