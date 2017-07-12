program posterfig3
use precision
implicit none
integer nmax,npt1,npt2,npt3
real, allocatable, dimension(:) :: rbb,px,py
real(double) :: temp
real(double), allocatable, dimension(:) :: time1,flux1,time2,flux2,     &
   time3,flux3,bb
character(80) :: file1,file2,file3

file1="phot1.dat"
file2="phot2.dat"
file3="phot3.dat"

nmax=1000
allocate(flux1(nmax),time1(nmax),flux2(nmax),time2(nmax),flux3(nmax),time3(nmax))
call readmod(file1,nmax,npt1,time1,flux1)
call readmod(file2,nmax,npt2,time2,flux2)
call readmod(file3,nmax,npt3,time3,flux3)
write(0,*) "npt: ",npt1,npt2,npt3
temp=flux1(1)
flux1=flux1/temp
temp=flux2(1)
flux2=flux2/temp
temp=flux3(1)
flux3=flux3/temp

!open PGPLOT device
call pgopen('?')!('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,2)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(2.5) !bigger text

allocate(bb(4))
bb(1)=minval(time1(1:npt1))
bb(2)=maxval(time1(1:npt1))
bb(3)=minval(flux1(1:npt1))
bb(4)=maxval(flux1(1:npt1))

allocate(rbb(4))
rbb=real(bb)
rbb(1)=0.0
rbb(2)=4.0
rbb(3)=rbb(3)-0.1*(rbb(4)-rbb(3))
rbb(4)=rbb(4)+0.1*(rbb(4)-rbb(3))

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4)) !plot scale
call pgbox("BCTS",0.0,0,"BCNTSV",0.0,0)
call pglabel("Time (hours)","","")
call pgptxt(rbb(1)-0.16*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2,90.0,0.5,          &
   "Relative Flux")

allocate(px(npt1),py(npt1))
px(1:npt1)=real(time1(1:npt1))
py(1:npt1)=real(flux1(1:npt1))
call pgline(npt1,px,py)
deallocate(px,py)

!allocate(px(npt2),py(npt2))
!px(1:npt2)=real(time2(1:npt2))
!py(1:npt2)=real(flux2(1:npt2))
!call pgsci(2)
!call pgline(npt2,px,py)
!call pgsci(1)
!deallocate(px,py)

allocate(px(npt3),py(npt3))
px(1:npt3)=real(time3(1:npt3))
py(1:npt3)=real(flux3(1:npt3))
call pgsci(2)
call pgline(npt3,px,py)
call pgsci(1)
deallocate(px,py)



end program posterfig3

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readmod(filename,nmax,npt,time,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nmax,npt
real(double), dimension(nmax) :: time,flux
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
   read(nunit,*,iostat=filestatus) flux(i)
   time(i)=dble(i)*2.5d0/60.0d0
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
