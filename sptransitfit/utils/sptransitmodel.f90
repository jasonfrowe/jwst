subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,     &
 exptime,ntt,tobs,omc,sptmodel)
use precision
implicit none
!import vars
integer :: nplanet,npars,nwv,nobs
integer, dimension(:) :: ntt
real(double), dimension(:) :: sol,time,exptime
real(double), dimension(:,:) :: solrange,sptmodel,tobs,omc
!local vars
integer :: nintg,iwv,ii,i,j
real(double), allocatable, dimension(:) :: tflux,bt
real(double) :: Pi,tPi,pid2,G,Cs,fDB,c1,c2,c3,c4,dil,voff,zpt,rhostar,  &
 epoch,per,b,rprs,ecw,esw,K,ted,ell,ag,bs2,eccn,w,adrs,incl,dnintg,     &
 tdnintg,dnintgm1,Eanom,phi0,ttcor,jm1,t,phi,Manom,Tanom,drs,x2,y2,     &
 trueanomaly,distance

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

!Model parameters
nintg=41 !number of samples to convolve integration time

!precompute doubles and repeatitive math
dnintg=dble(nintg) !convert integer to double
tdnintg=2.0d0*dnintg
dnintgm1=2.0d0*dnintg-2.0d0

!Physical Constants
Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi
pid2=Pi/2.0d0  !and Pi/2
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
Cs=2.99792458d8 !Speed of light
fDB=1.0d0 !Doppler Boosting factor

!init sptmodel to zero
sptmodel=0.0d0

!allocate array for calculating flux model and impact parameter over
!exposure time
allocate(tflux(nintg),bt(nintg))

do iwv=1,nwv !loop over all bandpasses

   !get parameters that do not depend on planet
   call getbasicpars(iwv,sol,solrange,rhostar,c1,c2,c3,c4,dil,voff,zpt)

   do ii=1,nplanet !loop over all planets

      call getplanetpars(iwv,ii,sol,solrange,epoch,per,b,rprs,ecw,esw,K,&
         ted,ell,ag)
      bs2=b*b !b^2
      call geteccn(ecw,esw,eccn,w) !get e,w given ecosw,esinw

      !calculate scaled semi-major axis from mean stellar density
      adrs=1000.0*rhostar*G*(Per*86400.0d0)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0)

      !mean,eccentric anomalies and phase offset for transit center
      Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
      Eanom=2.0d0*atan(Eanom)
      phi0=Eanom-eccn*sin(Eanom)

      !look over all time-steps for a single bandpass

      do i=1,nobs

         !transit-time variations
         call lininterp(tobs,omc,nplanet,nobs,ii,ntt,time(i),ttcor)

         !integrate over exposure time
         do j=1,nintg
            jm1=dble(j-1) !pre-compute double
            tflux(j)=0.0 !initialize model flux to zero

            !sample small dt across exposure time for integration.
            t=time(i)-exptime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg)-      &
               epoch-ttcor

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
            incl=acos(b/drs)
            x2=drs*Sin(Tanom-w)
            y2=drs*Cos(Tanom-w)*cos(incl)

            !time specific impact parameter.
            bt(j)=sqrt(x2*x2+y2*y2)

         enddo

         read(5,*) !simple pause statement

      enddo

   enddo

enddo


return
end subroutine sptransitmodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine geteccn(ecw,esw,eccn,w)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
real(double) :: ecw,esw,eccn,w
!local vars
real(double) :: Pi,tPi,pid2

Pi=acos(-1.d0) !define Pi
tPi=2.0d0*Pi   !and 2*Pi
pid2=Pi/2.0d0  !and Pi/2

eccn=(ecw*ecw+esw*esw)
if(eccn.ge.1.0) eccn=0.99
if(eccn.eq.0.0d0)then
   w=0.0d0
else
   if(ecw.eq.0.0d0)then
      w=pid2
   else
      w=atan(esw/ecw)
   endif
   if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
      w=tPi+w
   elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then
      w=Pi+w
   elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
      w=Pi+w
   endif
endif

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getplanetpars(iwv,nplanet,sol,solrange,epoch,per,b,rprs,ecw, &
 esw,K,ted,ell,ag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer iwv,nplanet
real(double) :: epoch,per,b,rprs,ecw,esw,K,ted,ell,ag
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solrange
!local vars
integer i,ii,j
integer, parameter :: npars=10
real(double), dimension(npars) :: lpars

j=0 !counter for indexing lpars
do i=9,8+npars
   j=j+1
   ii=i+10*(nplanet-1) !grab appropriate index for multi-planets
   if(solrange(ii,2)-solrange(ii,1).gt.0)then
      lpars(j)=sol(solrange(ii,1)+iwv-1)
   else
      lpars(j)=sol(solrange(ii,1))
   endif
enddo

epoch=lpars(1)
per=lpars(2)
b=abs(lpars(3)) !b is always positive
rprs=abs(lpars(4)) !Rp/R* is alway positive
ecw=lpars(6) !note the reversal here.
esw=lpars(5)
K=lpars(7)
ted=lpars(8)
ell=lpars(9)
ag=lpars(10)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getbasicpars(iwv,sol,solrange,rhostar,c1,c2,c3,c4,dil,voff,  &
 zpt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: iwv
real(double) :: rhostar,c1,c2,c3,c4,dil,voff,zpt
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solrange
!local vars
integer i
integer, parameter :: npars=8
real(double), dimension(npars) :: lpars

do i=1,npars
   if(solrange(i,2)-solrange(i,1).gt.0)then
      lpars(i)=sol(solrange(i,1)+iwv-1)
   else
      lpars(i)=sol(solrange(i,1))
   endif
enddo

rhostar=lpars(1) !mean stellar density
c1=lpars(2)      !limb-darkening
c2=lpars(3)
c3=lpars(4)
c4=lpars(5)
dil=lpars(6)     !dilution (0 == no dilution)
voff=lpars(7)    !velocity offset
zpt=lpars(8)     !photometric offset

return
end
