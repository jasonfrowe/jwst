!simple program to plot delta's between traces
program plotdtraces
use precision
implicit none
integer iargc,nunit,filestatus,i,j,k,nmax,npt,ntrace,ip,nfit,nplot
real, allocatable, dimension(:) :: px,py
real(double) :: dumr,dp,pwcpos(11)
real(double), allocatable, dimension(:) :: rpos,sol,pos
character(80) :: filename,dumc

pwcpos(1)=245.8617706298828
pwcpos(2)=245.8642272949219
pwcpos(3)=245.8568878173828
pwcpos(4)=245.8471069335938
pwcpos(5)=245.8495483398438
pwcpos(6)=245.8568878173828
pwcpos(7)=245.8544464111328
pwcpos(8)=245.8544464111328
pwcpos(9)=245.8666687011719
pwcpos(10)=245.8642272949219
pwcpos(11)=245.8495483398438

if(iargc().lt.2)then
   write(0,*) "Usage: trace1, trace2... "
   write(0,*) "  must provide at least 2 traces"
   stop
endif

call getarg(1,filename)
nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

read(nunit,*) ntrace,nmax !number of traces and number of entries
nfit=1+ntrace*9
nmax=nmax+1

!allocate space for reference pos
allocate(rpos(nmax),sol(nfit),pos(nmax),px(nmax),py(nmax))

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax,i
      stop
   endif
   read(nunit,*,iostat=filestatus) ip,(sol(j),j=1,10)
   if(filestatus == 0) then
      rpos(i)=sol(9)+(sol(3)+sol(6))/2.0d0
!      write(0,*) rpos(i),sol(9),sol(3),sol(6)
!      read(5,*)
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
close(nunit) !close file
npt=i-1
write(0,*) "npt: ",npt

call pgopen('?')  !'?' lets the user choose the device.
call PGPAP (8.0,1.0) !use a square 8" across
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgwindow(0.0,2048.0,-1.0,1.0)
CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
call pglabel("X (pixels)","dY (pixels)","")

do k=2,iargc()

   call getarg(k,filename)
   write(0,'(A80)') filename
   nunit=10
   open(unit=nunit,file=filename,iostat=filestatus,status='old')
   if(filestatus>0)then !trap missing file errors
      write(0,*) "Cannot open ",filename
      stop
   endif

   read(nunit,*) dumc !don't need header info anymore (could do checks)

   i=1
   do
      if(i.gt.nmax)then
         write(0,*) "Increase nmax to match data points"
         write(0,*) "nmax: ",nmax,i
         stop
      endif
      read(nunit,*,iostat=filestatus) ip,(sol(j),j=1,10)
      if(filestatus == 0) then
         pos(i)=sol(9)+(sol(3)+sol(6))/2.0d0
         i=i+1
      elseif(filestatus == -1) then
         exit  !successively break from data read loop.
      else
         write(0,*) "File Error!! Line:",i+1
         write(0,900) "iostat: ",filestatus
         stop
      endif
   enddo
   close(nunit) !close file
   npt=i-1
   write(0,*) "npt: ",npt

   nplot=0
   do j=1,npt
!      write(0,*) pos(j),rpos(j)
!      read(5,*)
      if((pos(j).gt.0.0).and.(rpos(j).gt.0.0))then
         dp=pos(j)-rpos(j)
         nplot=nplot+1
         px(nplot)=real(j)
         py(nplot)=real(dp)
      endif
   enddo
   call pgsci(k-1)
   call pgline(nplot,px,py)
enddo

call pgclos()

end program plotdtraces
