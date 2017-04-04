program limbdarkpars
use precision
use fittingmod
implicit none
integer nfit,info,lwa,nmu,iargc,filestatus,nunit,idum,i,nplot,nstep,    &
   nbin,j
integer, allocatable, dimension(:) :: iwa
real(double) :: tol,rwv,bandmu,banddmu,wv,dm1,dm2,norm
real(double), allocatable, dimension(:),target :: mu,dIn,rdIn
real(double), allocatable, dimension(:) :: sol,fvec,wa
real, allocatable, dimension(:) :: px,py
character(80) :: filename,cline,title
external fcn !source for fcn can be found at bottom of this file

nmu=17 !number of surface angles

!parameters for fitting limb-darkening..
nfit=2

!if nplot=0 then plot, else do not plot
nplot=0

allocate(sol(nfit),fvec(nmu),iwa(nfit))
if(nfit.eq.4)then
   sol(1)=0.0d0
   sol(2)=0.0d0
   sol(3)=0.0d0
   sol(4)=1.0d0
else !default it to use quadratic law
   nfit=2
   sol(1)=0.0d0
   sol(2)=1.0d0
endif

tol=1.0d-8
info=0
lwa=nmu*nfit+5*nmu*nfit
allocate(wa(lwa))

if(iargc().lt.3)then
   write(0,*) "Usage: limbdarkpars r3500_I.dat mu dmu"
   write(0,*) " r3500_I.dat - output from ATLAS-9 models"
   write(0,*) " mu - central wavelength for bandpass (um)"
   write(0,*) " dmu - +/- width of bandpass (um)"
   stop
endif

call getarg(1,filename)  !get filename for photometry
nunit=10 !unit number for file list
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

!read in location of center of bandpass
call getarg(2,cline)
read(cline,*) bandmu

!read in +/- width of bandpass
call getarg(3,cline)
read(cline,*) banddmu
if(banddmu.le.0.0d0)then
   write(0,*) "Error: dmu must be greater than zero"
   stop
endif

!assign wavelength range
dm1=(bandmu-banddmu)*10000.0d0 !um -> A
dm2=(bandmu+banddmu)*10000.0d0

!read in viewing angles
allocate(mu(nmu))
read(nunit,*) idum,(mu(i),i=1,nmu)

if(nplot.eq.0)then
   call pgopen('/xserve') !open PGPlot device
   call pgask(.false.) !false=don't ask for new page.. just do it.
   call PGPAP ( 6.0 ,1.0) !paper size
   call pgsubp(1,1)  !break up plot into grid
   call pgpage()
   call pgvport(0.20,0.9,0.20,0.9)
   call pgwindow(-0.1,1.1,-0.1,1.1)
   call pgsch(2.0)
endif

!read in wavelengths and intensities
allocate(dIn(nmu)) !allocate array for intensities
dIn=0.0d0
dIn2 => dIn !update pointer for fitting routine

if(nplot.eq.0) then
   allocate(px(nmu),py(nmu)) !allocate arrays for plotting
   do i=1,nmu
      px(i)=real(mu(i))
   enddo
endif

mu2 => mu  !update pointer for fitting routine

allocate(rdIn(nmu))
rwv=0.0d0
rdIn=0.0d0

nstep=0 !count number of steps
nbin=0  !number of flux measurements inside bandpass
do
   read(nunit,*,iostat=filestatus) wv,(dIn(i),i=1,nmu)
   wv=wv*10.0 !nm -> A
   dIn=4*dIn*2.99792458d18/(wv*wv) !convert to ergs s-1 cm-2 A-1
   if(filestatus == 0) then
      nstep=nstep+1 !count number of steps
      if((wv.ge.dm1).and.(wv.le.dm2))then
         rwv=rwv+wv    !sum up wavelength and flux to get averages
         rdIn=rdIn+dIn
         nbin=nbin+1
      endif
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",nstep+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

!calculate bandpass average
wv=rwv/dble(nbin)
dIn=rdIn/dble(nbin)

norm=dIn(1)
do i=1,nmu
   dIn(i)=dIn(i)/norm !normalize to center of star
   if(nplot.eq.0)then
      py(i)=real(dIn(i))
   endif
enddo

!plot the binned model fluxes
if(nplot.eq.0)then
   call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
   call pglabel("mu","Intensity","")
   call pgline(nmu,px,py)
   call pgpt(nmu,px,py,17)
   write(title,503) wv
   503 format(F9.3)
   call pgsch(1.5)
   call pgptxt(0.0,1.0,0.0,0.0,title)
   call pgsch(2.0)
endif

!initialize fit parameters
if(nfit.eq.4)then
   sol(1)=0.0d0
   sol(2)=0.0d0
   sol(3)=0.0d0
   sol(4)=1.0d0
else !default it to use quadratic law
   nfit=2
   sol(1)=0.24d0
   sol(2)=0.44d0
endif

!fit the data.  fcn subroutine can be found at end of this source
call lmdif1(fcn,nmu,nfit,sol,fvec,tol,info,iwa,wa,lwa)
!write(0,*) "info: ",info,wv
!write(6,500) wv,(sol(i),i=1,nfit),norm
write(6,500) (sol(i),i=1,nfit)
500 format(40(1X,1PE17.10))

if(nfit.eq.2)then
   do i=1,nmu
      py(i)=real(1.0d0-sol(1)*(1.0d0-mu(i))-sol(2)*(1.0d0-mu(i))**2.0d0)
   enddo
else
   do i=1,nmu
      py(i)=1.0
      do j=1,nfit
         py(i)=py(i)-real(sol(j)*(1.0-mu(i)**(dble(j)/2.0d0)))
      enddo
   enddo
endif
call pgsci(2)
call pgslw(3)
call pgline(nmu,px,py)
call pgslw(1)
call pgsci(1)

end program limbdarkpars

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(m,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: m,n,iflag,i,j
real(double), dimension(n) :: x
real(double), dimension(m) :: fvec

if(n.eq.2)then
   do i=1,m
      fvec(i)=1.0d0-x(1)*(1.0d0-mu2(i))-x(2)*(1.0d0-mu2(i))**2.0d0
   enddo
else
   do i=1,m
      fvec(i)=1.0d0
      do j=1,n
         fvec(i)=fvec(i)-x(j)*(1.0d0-mu2(i)**(dble(j)/2.0d0))
      enddo
   enddo
endif

do i=1,m
   fvec(i)=(fvec(i)-dIn2(i))/1.0d0
enddo

!priors to get physical parameters (quadratic only for the moment)

if(n.eq.2)then
   if((x(1)+x(2).gt.1.0d0).or.(x(1).lt.0.0d0).or.                       &
    (x(1)+2.0d0*x(2).lt.0.0))then
      fvec=99.9e30
   endif
endif


return
end
