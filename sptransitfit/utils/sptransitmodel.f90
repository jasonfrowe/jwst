subroutine sptransitmodel(nplanet,npars,sol,solrange,nwv,nobs,time,     &
 exptime,ntt,tobs,omc,sptmodel,nwvc)
use precision
implicit none
!import vars
integer :: nplanet,npars,nwv,nobs,nwvc
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: sptmodel,tobs,omc,time,exptime
!local vars
integer :: nintg,iwv,ii,i,j,caltran,nwv1,nwv2
real :: twork
real(double), allocatable, dimension(:) :: tflux,bt,vt,tide,alb,mu,bp,  &
 occ,time1
real(double) :: Pi,tPi,pid2,G,Cs,fDB,c1,c2,c3,c4,dil,voff,zpt,rhostar,  &
 epoch,per,b,rprs,ecw,esw,K,ted,ell,ag,eccn,w,adrs,incl,dnintg,         &
 tdnintg,dnintgm1,Eanom,phi0,ttcor,jm1,t,phi,Manom,Tanom,drs,x2,y2,     &
 trueanomaly,distance,albedomod,tm,ratio

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

!Model parameters
nintg=11 !number of samples to convolve integration time

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

!allocate array for calculating flux model, impact parameter, radial
!velocity, ellipodial variations and reflective/emission phase changes
!over exposure time
allocate(tflux(nintg),bt(nintg),vt(nintg),tide(nintg),alb(nintg),       &
 mu(nintg),bp(nintg),occ(nintg))
allocate(time1(nobs)) !needed for lininterp.

!added ability to calcalate likelihood for only a single bandpass
!this can significantly speed up calculations
if(nwvc.eq.0)then !loop over all bandpasses
   nwv1=1
   nwv2=nwv
else              !loop over single bandpass
   nwv1=nwvc
   nwv2=nwvc
endif

do iwv=nwv1,nwv2 !loop over all bandpasses

   !get parameters that do not depend on planet
   call getbasicpars(iwv,sol,solrange,rhostar,c1,c2,c3,c4,dil,voff,zpt)

   do ii=1,nplanet !loop over all planets

      call getplanetpars(iwv,ii,sol,solrange,epoch,per,b,rprs,ecw,esw,K,&
         ted,ell,ag)
      call geteccn(ecw,esw,eccn,w) !get e,w given ecosw,esinw

      !calculate scaled semi-major axis from mean stellar density
      adrs=1000.0*rhostar*G*(Per*86400.0d0)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0)

      !inclination
      incl=acos(b/adrs)

      !mean,eccentric anomalies and phase offset for transit center
      Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
      Eanom=2.0d0*atan(Eanom)
      phi0=Eanom-eccn*sin(Eanom)

      !look over all time-steps for a single bandpass

      time1(:)=time(iwv,:)
      do i=1,nobs

         !transit-time variations
         call lininterp(tobs,omc,nplanet,nobs,ii,ntt,time1(i),ttcor)

         !integrate over exposure time
         do j=1,nintg
            jm1=dble(j-1) !pre-compute double
            tflux(j)=0.0 !initialize model flux to zero

            !sample small dt across exposure time for integration.
            t=time(iwv,i)-exptime(iwv,i)*                               &
             (0.5d0-1.0d0/tdnintg-jm1/dnintg)-epoch-ttcor

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

            !time specific impact parameter.
            bt(j)=sqrt(x2*x2+y2*y2)

            !Radial velocity - used for Doppler Beaming
            vt(j)=K*(cos(Tanom-w+pid2)+eccn*cos(-w+pid2))

            !Ellipsodial variations
            tide(j)=ell*(drs/adrs)**(1.0d0/3.0d0)*                      &
             cos(2.0d0*(Pid2+Tanom-w))

            !Albedo/thermal phase changes
            alb(j)=albedomod(Pi,ag,Tanom-w)*adrs/drs

         enddo

         if(y2.ge.0.0d0)then  !If we have a potential transit
            caltran=0 !if zero, there is no transit
            !scan though all calculating bt and see if a transit occurs
            do j=1,nintg
               if(bt(j).le.1.0d0+RpRs)then !condition for a transit
                  caltran=1
               endif
            enddo
            if(caltran.eq.1) then
               !quadratic co-efficients
               if((c3.eq.0.0).and.(c4.eq.0.0))then
                  call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
               !Kipping co-efficients
               elseif((c1.eq.0.0).and.(c2.eq.0.0))then
                  c1=2.0d0*sqrt(c3)*c4 !convert to regular LD
                  c2=sqrt(c3)*(1.0d0-2.0d0*c4)
                  call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
                  c1=0.0d0  !zero out entries.
                  c2=0.0d0
               !non-linear law.
               else
                  !non-linear code is buggy.  Use with caution
                  !call occultnl(RpRs,c1,c2,c3,c4,bt,tflux,mulimbf,nintg)
                  !using small planet approx. for now.  May not be great
                  !for hot-Jupiters
                  call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
               endif

            else !there is no transit. So init flux to 1.0
               do j=1,nintg
                  tflux(j)=1.0d0
               enddo
            endif

            !now we sum up all the model fluxes across the exposure time
            tm=0.0d0
            do j=1,nintg
               if(RpRs.le.0.0) tflux(j)=1.0d0 !make sure r/R* is sane
!                 model=transit+doppler+ellipsodial
               tm=tm+tflux(j)-fDB*vt(j)/Cs+tide(j)+alb(j)
            enddo
            tm=tm/dnintg !calculate average flux across integration time

         else !We have an eclipse/occultation
            tm=0.0d0
            do j=1,nintg
               bp(j)=bt(j)/RpRs !impact parameter of eclipse/occultation
            enddo
            call occultuniform(bp,1.0/RpRs,occ,nintg)
            !integrate over exposure time
            do j=1,nintg
               ratio=1.0d0-occ(j) !scaling to modeled occultation depth
               if(RpRs.le.0.0d0) ratio=0.0d0
               !add in occultation, doppler, tidal and phase curve
               tm=tm+(1.0d0-ted*ratio)-fDB*vt(j)/Cs+tide(j)+alb(j)
            enddo
            tm=tm/dnintg !average model flux over exposure time
         endif !end of if statement for transit vs occult.

         tm=tm+(1.0d0-tm)*dil !add dilution

         !add flux for this planet to the overall model
         sptmodel(iwv,i)=sptmodel(iwv,i)+tm
         !write(6,'(2050(1PE17.10,1X))') time1(i),tm,(bt(j),j=1,nintg)

      enddo !end of looping over all time steps
      !read(5,*)

   enddo !end of looping over all planets

   !adding in zeropoint.
   sptmodel(iwv,:)=zpt*sptmodel(iwv,:)

enddo !end of looping over all bandpasses


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
integer, dimension(:,:) :: solrange
real(double) :: epoch,per,b,rprs,ecw,esw,K,ted,ell,ag
real(double), dimension(:) :: sol
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
integer, dimension(:,:) :: solrange
real(double) :: rhostar,c1,c2,c3,c4,dil,voff,zpt
real(double), dimension(:) :: sol
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
