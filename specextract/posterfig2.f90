program posterfig2
use precision
implicit none
integer nmax,npt,nplot,npt2,npt3,npt4
real :: t1,t2
real, allocatable, dimension(:) :: rbb,px,py
real(double), allocatable, dimension(:) :: wv,flux,bb,wv2,flux2,wv3,flux3, &
   wv4,flux4
character(80) :: filename,filename2,filename3

nmax=15000
filename="spgen_c_b20.trace"
filename2="spgen_c_b00.trace"
filename3="ref.txt"

allocate(wv(nmax),flux(nmax))
call readflux(filename,nmax,npt,wv,flux)
write(0,*) "Number of data points read: ",npt

!open PGPLOT device
call pgopen('?')!('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,0.7) !use a square 8" across
call pgsubp(1,3)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(2.5) !bigger text

allocate(bb(4))
bb(1)=minval(wv(1:npt))
bb(2)=maxval(wv(1:npt))
bb(3)=minval(flux(1:npt))
bb(4)=maxval(flux(1:npt))

allocate(rbb(4))
rbb=real(bb)
rbb(4)=rbb(4)+0.1*(rbb(4)-rbb(3))

call pgvport(0.15,0.95,0.0,0.95) !make room around the edges for labels
call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4)) !plot scale
call pgbox("BCTS",0.0,0,"BCNTSV",0.0,0)
!call pglabel("","Flux (e-)","")
call pgptxt(rbb(1)-0.08*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2,90.0,0.5,          &
   "Flux (e-)")

allocate(px(npt),py(npt))
px(1:npt)=real(wv(1:npt))
py(1:npt)=real(flux(1:npt))

call pgline(npt,px,py)
deallocate(px,py)

allocate(wv2(nmax),flux2(nmax))
call readflux(filename2,nmax,npt2,wv2,flux2)
!flux2=flux2/1.01467
write(0,*) "Number of data points read: ",npt2

allocate(px(npt2),py(npt2))
px(1:npt2)=real(wv2(1:npt2))
py(1:npt2)=real(flux2(1:npt2))

call pgsci(2)
call pgline(npt2,px,py)
call pgsci(1)

call pgpage()
py(1:npt)=real(flux2(1:npt)/flux(1:npt))

rbb(3)=minval(py(1:npt))
rbb(4)=maxval(py(1:npt))
t1=rbb(3)
t2=rbb(4)
rbb(3)=0.975 !t1-0.1*(t2-t1)
rbb(4)=0.985 !t2+0.1*(t2-t1)

call pgvport(0.15,0.95,0.15,1.00) !make room around the edges for labels
call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4)) !plot scale
call pgbox("BCNTS",0.0,0,"BCNTSV",0.0,0)
call pglabel("Wavelength (um)","","")
call pgptxt(rbb(1)-0.08*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2,90.0,0.5,          &
   "Ratio (e-)")

call pgline(npt,px,py)

allocate(wv3(nmax),flux3(nmax))
call readref(filename3,nmax,npt3,wv3,flux3)
flux3=flux3*0.99865
write(0,*) "npt3: ",npt3

px(1:npt3)=real(wv3(1:npt3))
py(1:npt3)=real(flux3(1:npt3))
call pgsci(2)
call pgline(npt3,px,py)
call pgsci(1)

deallocate(px,py)
allocate(wv4(nmax),flux4(nmax))
call readmod("hd209458-rprs-fix.txt",nmax,npt4,wv4,flux4)
write(0,*) "npt4: ",npt4
allocate(px(npt4),py(npt4))
px(1:npt4)=real(wv4(1:npt4))
py(1:npt4)=real(flux4(1:npt4))
call pgsci(3)
call pgline(npt4,px,py)
call pgsci(1)


call pgclos()

end program posterfig2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readmod(filename,nmax,npt,wv,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nmax,npt
real(double), dimension(nmax) :: wv,flux
character(80) :: filename
!local vars
integer :: nunit,filestatus,i,dumc,j
real(double) :: dumr,p2w,wt

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) wv(i),flux(i)
!   wt=p2w(wv(i),1,1)
   wv(i)=wv(i)/10000.0
   flux(i)=1.0-flux(i)*flux(i)
!   write(0,*) wv(i),flux(i)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1


return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readref(filename,nmax,npt,wv,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nmax,npt
real(double), dimension(nmax) :: wv,flux
character(80) :: filename
!local vars
integer :: nunit,filestatus,i,dumc,j
real(double) :: dumr,p2w,wt

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) wv(i),flux(i)
!   wt=p2w(wv(i),1,1)
   wv(i)=wv(i)/10000.0
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1


return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readflux(filename,nmax,npt,wv,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nmax,npt
real(double), dimension(nmax) :: wv,flux
character(80) :: filename
!local vars
integer :: nunit,filestatus,i,dumc,j
real(double) :: dumr,p2w,wt

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

read(nunit,*,iostat=filestatus) dumc !first line is not needed.

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) wv(i),(dumr,j=1,7),flux(i)
   wt=p2w(wv(i),1,1)
   wv(i)=wt/10000.0
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1


return
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function p2w(p,noversample,ntrace)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!convert pixel to wavelength
! n=1 polynomial fit
use precision
implicit none
integer :: i,noversample,ntrace
integer, parameter :: nc=5
real(double) :: pix,p,p2w
real(double), dimension(nc) :: c1,c2,c3,c
data c1/2.60188,-0.000984839,3.09333e-08,-4.19166e-11,1.66371e-14/
data c2/1.30816,-0.000480837,-5.21539e-09,8.11258e-12,5.77072e-16/
data c3/0.880545,-0.000311876,8.17443e-11,0.0,0.0/

select case(ntrace)
   case(1)
      c=c1
   case(2)
      c=c2
   case(3)
      c=c3
end select

pix=p/dble(noversample)
p2w=c(1)
do i=2,nc
   p2w=p2w+pix**dble(i-1)*c(i)
enddo
p2w=p2w*10000.0 !um->A

!write(0,*) p2w,pix,p
!read(5,*)

return
end
