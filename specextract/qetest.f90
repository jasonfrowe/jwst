program qetest
use precision
implicit none
integer nmax,nunit,nf,filestatus,i,j,k
integer, allocatable, dimension(:) :: npt
real(double) :: w,f,y
real(double), allocatable, dimension(:) :: y2,y3,a
real(double), allocatable, dimension(:,:) :: wv,ff
character(80) :: dumc
character(80), allocatable, dimension(:) :: file

nmax=7000 !maximum number of lines in a file
nf=3 !number of files
!allocate space to read in wavelengths and "fluxes"
allocate(wv(nf+2,nmax),ff(nf+2,nmax))
!allocate integer to track number of data points read in
allocate(npt(nf+2))
!allocate space for filenames
allocate(file(nf))
file(1)="2928_TUNG.dc.ap" !extracter aperture from CV3 data Tungsten obs
file(2)="otp714_rehearsal_osim_only.dat" !reference Tungsten
file(3)="thr_gr700_o1_alone.txt" !expected blaze function

!co-efficients to convert from pixel to wavelength
allocate(a(3))
a(1)=-905.219d0
a(2)=1.09890d0
a(3)=-2.03489d-5

nunit=10

do i=1,nf
   open(unit=nunit,file=file(i),iostat=filestatus,status='old')
   if(filestatus>0)then !trap missing file errors
      write(0,*) "Cannot open ",file(i)
      stop
   endif

!  ignore first line of file 2.
   if(i.eq.2) read(nunit,*) dumc

   j=0
   do
      if(j.gt.nmax)then
         write(0,*) "Increase nmax to match data points"
         write(0,*) "nmax: ",nmax,i
         stop
      endif
      read(nunit,*,iostat=filestatus) w,f
      if(filestatus == 0) then
         if(f.ge.0.0d0)then
            j=j+1
            if(i.eq.2)then
               wv(i,j)=w*1000.0d0 !um to nm
            elseif(i.eq.1)then
               wv(i,j)=(sqrt(-4.0*a(1)*a(3)+a(2)*a(2)+4.0*a(3)*w)-a(2))/(2.0*a(3))
            else
               wv(i,j)=w
            endif
            ff(i,j)=f
         endif
!         write(0,*) j,wv(i,j),ff(i,j)
      elseif(filestatus == -1) then
         exit  !successively break from data read loop.
      else
         write(0,*) "File Error!! Line:",i+1
         write(0,900) "iostat: ",filestatus
         900 format(A8,I3)
!         stop
      endif
   enddo
   npt(i)=j !store number of data points read in
   write(0,*) "npt",i,npt(i)
   close(nunit)
enddo

!now we need to set up a spline for files 2 and 3
allocate(y2(npt(2)))
call spline(wv(2,1:npt(2)),ff(2,1:npt(2)),npt(2),1.d30,1.d30,y2)
allocate(y3(npt(3)))
call spline(wv(3,1:npt(3)),ff(3,1:npt(3)),npt(3),1.d30,1.d30,y3)

!now fill in values
do i=1,npt(1)
   call splint(wv(2,1:npt(2)),ff(2,1:npt(2)),y2,npt(2),wv(1,i),y)
   wv(4,i)=wv(1,i)
   ff(4,i)=y
   call splint(wv(3,1:npt(3)),ff(3,1:npt(3)),y3,npt(3),wv(1,i),y)
   wv(5,i)=wv(1,i)
   ff(5,i)=y
   f=ff(4,i)*ff(5,i)/(ff(2,i))
   write(6,'(5(1PE17.10,1X))') wv(1,i),ff(1,i),ff(4,i),ff(5,i),f
enddo

end program qetest
