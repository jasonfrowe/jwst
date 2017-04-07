subroutine fittransitmodel8(npars,sol,solerr,solrange,nwv,nobs,time,    &
 flux,ferr,exptime,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: npars,nwv,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc, &
 solrange
!local vars
integer :: n,i
real(double) :: tol
real(double), allocatable, dimension(:) :: solin

interface
   subroutine EstZpt(npars,sol,solerr,solrange,nwv,nobs,time,flux,      &
    exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: npars,nwv,nobs
      integer, dimension(:) :: ntt
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc,&
       solrange
   end subroutine EstZpt
end interface

!Get estimate for zero points.
call EstZpt(npars,sol,solerr,solrange,nwv,nobs,time,flux,exptime,ntt,   &
 tobs,omc)

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
