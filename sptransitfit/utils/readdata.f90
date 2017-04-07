!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readdata(nunit,nobsmax,nwvmax,nobs,nwv,time,flux,ferr,       &
 exptime)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!(c) Jason Rowe 2017
use precision
implicit none
!input vars
integer nunit,nobsmax,nwvmax,nobs,nwv
real(double), dimension(nwvmax,nobsmax) :: time,flux,ferr,exptime
!local vars
integer i,j,filestatus
real(double), allocatable, dimension(:) :: t,f,fe,et
character(80) dumc

!first line contains info about number of bandpasses that we already
!read in
read(nunit,*) dumc
!second line contains labels for the bandpasses.
read(nunit,*) dumc

!t and f are temp variables for reading in data and making sure we do
!not run into errors
allocate(t(nwvmax),f(nwvmax),fe(nwvmax),et(nwvmax))

i=0 !counter for number of observations.
do
   if(i.gt.nobsmax)then
      write(0,*) "nobsmax is too small. This should not happen!"
      write(0,*) "nobsmax: ",nobsmax,i
      stop
   endif
   read(nunit,*,iostat=filestatus) (t(j),f(j),fe(j),et(j),j=1,nwvmax)
   if(filestatus == 0) then
      i=i+1
      time(:,i)=t(:) !assign observation time to array
      flux(:,i)=f(:) !assign flux measurements to array
      ferr(:,i)=fe(:)
      exptime(:,i)=et(:)/86400.0d0 !convert seconds to days to match time(i)
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

nwv=nwvmax
nobs=i

return
end subroutine readdata

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getnobsnwv(nunit,nobsmax,nwvmax)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!(c) Jason Rowe 2017
use precision
implicit none
!input vars
integer nunit,nobsmax,nwvmax
!local vars
integer i,filestatus
real(double) dumr
character(80) dumc

!first line contains the number of wavelengths
read(nunit,*) dumc,nwvmax
!second line contains labels for the bandpasses
read(nunit,*) dumc

!rest of lines are the time and flux measurements.  We just need to
!count the number of lines.

i=0 !counter for number of lines
do
   read(nunit,*,iostat=filestatus) dumr !first column is obstime
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
nobsmax=i

!rewind file so we are at the top again
rewind(nunit)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine openfile(nunit,obsfile)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!(c) Jason Rowe 2017
!input vars
implicit none
integer nunit
character(80) obsfile
!local vars
integer filestatus

nunit=10
open(unit=nunit,file=obsfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",obsfile
   stop
endif

return
end

