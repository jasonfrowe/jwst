subroutine EstZpt(npars,nplanet,sol,solerr,solrange,nwv,nobs,time,flux, &
exptime,ntt,tobs,omc)
use precision
implicit none
integer :: npars,nwv,nobs,nplanet
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,exptime,tobs,omc,      &
   solrange
!local vars
real(double), allocatable, dimension(:,:) :: bt

!deal with case of no out-of-transit data.

interface
   subroutine getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,   &
    omc,bt)
      use precision
      implicit none
      integer :: nwv,nplanet,npars,nobs
      integer, dimension(:) :: ntt
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: bt,solrange,time,tobs,omc
   end subroutine getb
end interface

!calculate if data is in transit or out of transit.
!this means calculating b and using rprs.
allocate(bt(nwv,nobs))
write(0,*) "Calling getb"
call getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,omc,bt)

!check if we have one zpt per bandpass, or just one.
if(solrange(8,2)-solrange(8,1).gt.0)then !8 corresponds to zero-point



else !else calculate for each bandpass

endif

return
end subroutine estzpt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,omc,  &
   bt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nwv,nplanet,npars,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: bt,solrange,time,tobs,omc
!local vars
integer :: i,ii,iwv
real(double), allocatable, dimension(:) :: time1
real(double) :: rhostar,c1,c2,c3,c4,dil,voff,zpt,epoch,per,b,rprs,ecw,  &
 esw,K,ted,ell,ag,eccn,w,adrs,Pi,tPi,pid2,G,Eanom,phi0,ttcor,t,phi,     &
 Manom,Tanom,trueanomaly,drs,incl,x2,y2,distance

!deal with transit-timing variations
!deal with multiple-transiting planets - return minimum value of b.

interface
   subroutine getbasicpars(iwv,sol,solrange,rhostar,c1,c2,c3,c4,dil,    &
    voff,zpt)
      use precision
      implicit none
      integer :: iwv
      real(double) :: rhostar,c1,c2,c3,c4,dil,voff,zpt
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solrange
   end subroutine getbasicpars

   subroutine getplanetpars(iwv,nplanet,sol,solrange,epoch,per,b,rprs,  &
    ecw,esw,K,ted,ell,ag)
      use precision
      implicit none
      integer iwv,nplanet
      real(double) :: epoch,per,b,rprs,ecw,esw,K,ted,ell,ag
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solrange
   end subroutine getplanetpars
end interface

!Physical Constants
Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi
pid2=Pi/2.0d0  !and Pi/2
G=6.674d-11 !N m^2 kg^-2  Gravitation constant

!Get value of impact parameter (bt) for all obsevation times.
!initiate bt to a large value
bt=9.9d30

allocate(time1(nobs)) !needed for lininterp.

do iwv=1,nwv !loop over all bandpasses

   !get parameters that do not depend on planet
   call getbasicpars(iwv,sol,solrange,rhostar,c1,c2,c3,c4,dil,voff,zpt)

   do ii=1,nplanet !loop over all planets

      !get parameters specific to one planet
      call getplanetpars(iwv,ii,sol,solrange,epoch,per,b,rprs,ecw,esw,K,&
         ted,ell,ag)
      call geteccn(ecw,esw,eccn,w) !get e,w given ecosw,esinw

      !calculate scaled semi-major axis from mean stellar density
      adrs=1000.0*rhostar*G*(Per*86400.0d0)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0)

      !mean,eccentric anomalies and phase offset for transit center
      Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
      Eanom=2.0d0*atan(Eanom)
      phi0=Eanom-eccn*sin(Eanom)

      !look over all time-steps for a single bandpass

      time1(:)=time(iwv,:)
      do i=1,nobs

         !transit-time variations
         call lininterp(tobs,omc,nplanet,nobs,ii,ntt,time1(i),ttcor)

         t=time(iwv,i)-epoch-ttcor
         !get orbital position (mean anomaly)
         phi=t/per-floor(t/per)
         phi=phi*tPi+phi0
         Manom=phi
         if(Manom.gt.tPi) Manom=Manom-tPi
         if(Manom.lt.0.0d0) Manom=Manom+tPi
         call kepler(Manom,Eanom,eccn)
         Tanom=trueanomaly(eccn,Eanom)
         if(phi.gt.Pi) phi=phi-tPi
         drs=distance(adrs,eccn,Tanom)
         incl=acos(b/drs) !angle at specific time step. (not orbital incl)
         x2=drs*Sin(Tanom-w)
         y2=drs*Cos(Tanom-w)*cos(incl)

         bt(iwv,i)=min(bt(iwv,i),sqrt(x2*x2+y2*y2)) !report min value of bt.
         !write(0,*) "bt:",time(iwv,i),bt(iwv,i)
         !read(5,*)

      enddo !end loop over all observations

   enddo !end loop over all planets

enddo !end loop over all bandpasses

return
end subroutine getb
