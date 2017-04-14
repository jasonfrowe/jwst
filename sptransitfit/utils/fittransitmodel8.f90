subroutine fittransitmodel8(npars,nplanet,sol,solerr,solrange,nwv,nobs, &
 time,flux,ferr,exptime,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: npars,nwv,nobs,nplanet
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc
!local vars
integer :: n,m,i,j,iprint,isave(44),iter
integer, allocatable, dimension(:) :: nbd,iwa
real :: twork
real(double) :: tol,factr,pgtol,dsave(29),f
real(double), allocatable, dimension(:) :: solin,l,u,g,wa,sol1
logical :: lsave(4)
character(60) :: task,csave
!export fit
integer :: nunit,filestatus
character(80) :: newfitfile

interface
   subroutine EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,   &
    flux,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: npars,nwv,nobs,nplanet
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc
   end subroutine EstZpt
   subroutine setbounds(n,nbd,l,u,nplanet,npars,solerr,solrange)
      use precision
      implicit none
      integer :: n,npars,nplanet
      integer, dimension(:) :: nbd
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: l,u
      real(double), dimension(:,:) :: solerr
   end subroutine setbounds
   function loglikelihood(nwv,nobs,nplanet,npars,sol,solrange,time,     &
    flux,ferr,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: time,flux,ferr,exptime,tobs,omc
      real(double) :: loglikelihood
   end function loglikelihood
   subroutine gradient(nwv,nobs,nplanet,npars,sol,solerr,solrange,time, &
    flux,ferr,exptime,ntt,tobs,omc,f,g)
      use precision
      implicit none
      integer :: nwv,nobs,nplanet,npars
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double) :: f
      real(double), dimension(:) :: sol,g
      real(double), dimension(:,:) :: solerr,time,flux,ferr,   &
       exptime,tobs,omc
   end subroutine gradient
   subroutine exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
      use precision
      implicit none
      integer :: nunit,npars,nplanet
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr
   end subroutine exportfitpars
end interface

!Get estimate for zero points.
call EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,flux,       &
 exptime,ntt,tobs,omc)

!export current fit
newfitfile="newfit.dat"
nunit=10
open(unit=nunit,file=newfitfile,iostat=filestatus)
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",newfitfile
   stop
endif
call exportfitpars(nunit,npars,nplanet,sol,solerr,solrange)
close(nunit)

!pick off variables that are being fitted.
 !solin contains only model parameters that are fitted
 !sol1 contains full set of updated model parameters
allocate(solin(npars),sol1(npars))
n=0 !number of parameters that are fit
do i=1,npars
   if(solerr(i,1).ne.0.0)then
      n=n+1
      solin(n)=sol(i) !contains subset of sol that is passed to ldif
   endif
enddo

allocate(nbd(n),l(n),u(n)) !allocate arrays that set bounds for fitted parameters
nbd=0 !default is that parameters are unbounded.
!set bounds for parameters, source for subroutine is in this file.
call setbounds(n,nbd,l,u,nplanet,npars,solerr,solrange)

allocate(g(n)) !contains gradient information
factr=1.0d+7 !1.d+12 for low, 1.d+7 for moderate, 1.d+1 for high accuracy
pgtol=1.0d-5 !projected gradient tolerance
m=5  !maximum number of variable metric corrections
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) ) !working array
allocate ( iwa(3*n) )  !integer working array
iprint=-1 !frequency and type of output generated

!setting up loop conditions for lbfgsb
task = 'START'

!loop for lbfgsb
iter=1 !count number of iterations
do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
               task.eq.'START')

   call setulb ( n, m, solin, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )

   write(0,'(A6,A6)') "task: ",task

   if (task(1:2) .eq. 'FG') then
      j=0
      do i=1,npars
         if(solerr(i,1).ne.0.0)then
            j=j+1
            sol1(i)=solin(j) !if parameter if being fitted
         else
            sol1(i)=sol(i) !parameter is fixed, grab from master
         endif
      enddo

      !calculate log(likelihood) with current model solution
      write(6,*) "Calling loglikelihood"
      if(iter.eq.1)then
         CALL CPU_TIME(twork)
         write(0,*) "Fstart: ",twork
      endif
      f=-loglikelihood(nwv,nobs,nplanet,npars,sol1,solrange,time,flux,   &
       ferr,exptime,ntt,tobs,omc)
      CALL CPU_TIME(twork)
      write(0,*) "F: ",f,twork

      !calculate gradient
      write(0,*) "Calculating g"
      call gradient(nwv,nobs,nplanet,npars,sol1,solerr,solrange,time,    &
       flux,ferr,exptime,ntt,tobs,omc,f,g)
      write(0,*) "G1: ",g(1)


      newfitfile="newfit.dat"
      nunit=10
      open(unit=nunit,file=newfitfile,iostat=filestatus)
      if(filestatus>0)then !trap missing file errors
         write(0,*) "Cannot open ",newfitfile
         stop
      endif
      call exportfitpars(nunit,npars,nplanet,sol1,solerr,solrange)
      close(nunit)


   endif

   iter=iter+1 !increase counter for number of iterations
!   read(5,*)

enddo

sol=sol1 !update solution.

write(0,*) "Hey.. we made it!"

return
end subroutine fittransitmodel8

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine setbounds(n,nbd,l,u,nplanet,npars,solerr,solrange)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: n,npars,nplanet
integer, dimension(:) :: nbd
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: l,u
real(double), dimension(:,:) :: solerr
!local vars
integer :: i,j,k,ii,np

k=0
do i=1,8 !global parameters
   if(i.eq.1)then !bounds on mean-stellar density (0 < rhostar)
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=1 !lower bound only
            l(k)=0.0d0 !lower bound
         endif
      enddo
   elseif(i.eq.2)then ! limb-darkening #1  0 < NL1
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=1 !lower bound only
            l(k)=0.0d0 !lower bound
         endif
      enddo
   elseif(i.eq.3)then ! limb-darkening #2 - unbounded (handled by prior)
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=0 !lower bound only
         endif
      enddo
   elseif(i.eq.4)then ! limb-darkening #3  0 < NL3 < 1
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=2 !lower bound only
            l(k)=0.0d0 !lower bound
            u(k)=1.0d0 !upper bound
         endif
      enddo
   elseif(i.eq.5)then ! limb-darkening #3  0 < NL4 < 1
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=2 !lower bound only
            l(k)=0.0d0 !lower bound
            u(k)=1.0d0 !upper bound
         endif
      enddo
   elseif(i.eq.6)then !dilution (0 < DIL < 1)
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=2 !both upper and lower bound
            l(k)=0.0d0 !lower bound
            u(k)=1.0d0 !upper bound
         endif
      enddo
   elseif(i.eq.7)then !velocity offset VOF (unbounded)
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=0 !unbounded.
         endif
      enddo
   elseif(i.eq.8)then !photometric offset ZPT (unbounded)
      do j=solrange(i,1),solrange(i,2)
         if(solerr(j,1).ne.0.0)then
            k=k+1
            nbd(k)=0 !unbounded.
         endif
      enddo
   endif
enddo

!planet parameters
do np=1,nplanet !loop over each planet in transit model
   do i=1,10    !loop over each parameter for a single planet
     ii=8+i+nplanet*(np-1)
      if(i.eq.1)then !EPO (unbounded)
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      elseif(i.eq.2)then !PER (unbounded)
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      elseif(i.eq.3)then !BB (0 < BB)
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=1 !lower bound only
               l(k)=0.0d0  !lower bound
            endif
         enddo
      elseif(i.eq.4)then !RDR (0 < RDR)
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=1 !lower bound only
               l(k)=0.0d0  !lower bound
            endif
         enddo
      elseif((i.eq.5).or.(i.eq.6))then !ECW/ESW (0 < ESW < 1)
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=2 !lower and upper bound
               l(k)=0.0d0  !lower bound
               u(k)=1.0d0  !upper bound - there is also an e<1 bound in likelihood
            endif
         enddo
      elseif(i.eq.7)then !KRV (unbounded) !radial velocity amplitude
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      elseif(i.eq.8)then !TED (unbounded) !occultation depth
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      elseif(i.eq.9)then !ELL (unbounded) !amplitude of ellipsodial variations
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      elseif(i.eq.10)then !ALB (unbounded) - amplitude of phase curve
         do j=solrange(ii,1),solrange(ii,2)
            if(solerr(j,1).ne.0.0)then
               k=k+1
               nbd(k)=0 !unbounded
            endif
         enddo
      endif
   enddo
enddo

write(0,*) 'k: ',k

return
end
