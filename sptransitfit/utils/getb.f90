!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getb(nwv,nplanet,npars,sol,solrange,nobs,time,ntt,tobs,omc,  &
   bt,imarktrans)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!get impact parameter for each observation, bandpass and planet (min)
!marks in-transit (imarktrans=0) and out of transit data (imarktrans=1)
use precision
implicit none
!import vars
integer :: nwv,nplanet,npars,nobs
integer, dimension(:) :: ntt
integer, dimension(:,:) :: imarktrans,solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: bt,time,tobs,omc
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
      integer, dimension(:,:) :: solrange
      real(double) :: rhostar,c1,c2,c3,c4,dil,voff,zpt
      real(double), dimension(:) :: sol
   end subroutine getbasicpars
   subroutine getplanetpars(iwv,nplanet,sol,solrange,epoch,per,b,rprs,  &
    ecw,esw,K,ted,ell,ag)
      use precision
      implicit none
      integer iwv,nplanet
      integer, dimension(:,:) :: solrange
      real(double) :: epoch,per,b,rprs,ecw,esw,K,ted,ell,ag
      real(double), dimension(:) :: sol
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

         if(bt(iwv,i).le.1.0d0+RpRs)then
            imarktrans(iwv,i)=0 !in-transit data
         else
            imarktrans(iwv,i)=1 !out-of-transit data
         endif

      enddo !end loop over all observations

   enddo !end loop over all planets

enddo !end loop over all bandpasses

return
end subroutine getb
