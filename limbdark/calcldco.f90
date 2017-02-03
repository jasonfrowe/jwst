program calcldco
!calculate limb-darkening coefficients
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
use fittingmod
implicit none
integer :: iargc,nunit,filestatus,nmu,idum,i,nstep,nfit,info,lwa,j,nskip,nbin
integer, allocatable, dimension(:) :: iwa
real(double) :: wv,norm,tol,dnbin,rwv
real(double), allocatable, dimension(:),target :: mu,dIn,rdIn
real(double), allocatable, dimension(:) :: sol,fvec,wa
real, allocatable, dimension(:) :: px,py
character(18) :: pname
character(80) :: filename,title
external fcn

nskip=1000 !frames to skip for plotting
nbin=1
dnbin=dble(nbin)

nmu=17 !number of surface angles

!parameters for fitting limb-darkening..
nfit=4
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

if(iargc().lt.1)then
   write(0,*) "Usage: calcldco r3500_I.dat"
   stop
endif

call getarg(1,filename)  !get filename for photometry
nunit=10 !unit number for file list
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

!read in viewing angles
allocate(mu(nmu))
read(nunit,*) idum,(mu(i),i=1,nmu)

call pgopen('/xserve') !open PGPlot device
call pgask(.false.) !false=don't ask for new page.. just do it.
call PGPAP ( 6.0 ,1.0) !paper size
call pgsubp(1,1)  !break up plot into grid
call pgpage()
call pgvport(0.20,0.9,0.20,0.9)
call pgwindow(-0.1,1.1,-0.1,1.1)
call pgsch(2.0)

!read in wavelengths and intensities
allocate(dIn(nmu)) !allocate array for intensities
dIn=0.0d0
dIn2 => dIn !update pointer for fitting routine
allocate(px(nmu),py(nmu)) !allocate arrays for plotting
do i=1,nmu
   px(i)=real(mu(i))
enddo

mu2 => mu  !update pointer for fitting routine

allocate(rdIn(nmu))
rwv=0.0d0
rdIn=0.0d0

nstep=0 !count number of steps
do
   read(nunit,*,iostat=filestatus) wv,(dIn(i),i=1,nmu)
   wv=wv*10.0 !nm -> A
   dIn=4*dIn*2.99792458d18/(wv*wv) !convert to ergs s-1 cm-2 A-1
   if(filestatus == 0) then
      nstep=nstep+1 !count number of steps
      write(pname,'(A1,I9.9,A8)') "p",nstep,".png/png"
      if(mod(nstep,nbin).ne.0)then
         rwv=rwv+wv
         rdIn=rdIn+dIn
         cycle
      else
         wv=(rwv+wv)/dnbin
         dIn=(rdIn+dIn)/dnbin
         rwv=0.0d0
         rdIn=0.0d0
      endif
      norm=dIn(1)
      do i=1,nmu
         dIn(i)=dIn(i)/norm !normalize to center of star
         py(i)=real(dIn(i))
!         write(6,*) py(i)
      enddo
      !plot the data
!      call pgpage()

      if(mod(nstep,nskip).eq.0)then
!         call pgopen(pname)
!         call PGPAP ( 8.0 ,1.0) !paper size
!         call pgsubp(1,1)  !break up plot into grid
!         call pgvport(0.20,0.9,0.20,0.9)
!         call pgwindow(-0.1,1.1,-0.1,1.1)
!         call pgsch(2.0)

         call pgeras()
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


      !fit the data
      call lmdif1(fcn,nmu,nfit,sol,fvec,tol,info,iwa,wa,lwa)
!      write(0,*) "info: ",info,wv
      write(6,500) wv,(sol(i),i=1,nfit),norm
      500 format(40(1X,1PE17.10))
      !plot the fit

      if(mod(nstep,nskip).eq.0)then
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
!      read(5,*)
!         call pgeras()
!         call sleepqq(30)
!         call pgclos()
      endif
      cycle
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",nstep+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
!call pgclos()

end program calcldco

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

return
end


