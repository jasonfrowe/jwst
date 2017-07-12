subroutine mhgmcmc(nplanet,npars,sol,solrange,nwv,nobs,time,exptime,ntt,tobs, &
 omc,seed)
!M-H algorithm with a Gibb's sampler 
use precision
implicit none
!import vars
integer :: nplanet,npars,nwv,nobs,seed
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: tobs,omc,time,exptime
!local vars
integer :: nwvc,n
real(double) :: ran2
real(double), allocatable, dimension(:) :: sol2
real(double), allocatable, dimension(:,:) :: sptmodel

allocate(sol2(npars)) !allocate array to make copy of model solution
sol2(1:npars)=sol(1:npars) !make copy of solution array

n=int(ran2(seed)*npars+1) !select a random parameter to vary.  

return
end subroutine mhgmcmc