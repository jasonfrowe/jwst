subroutine plotrdr(nwv,rdr)
use precision
implicit none
!import vars
integer :: nwv
real(double), dimension(:) :: rdr
!local vars
real :: bb(4),i
real, allocatable, dimension(:) :: xp,yp
real(double) :: minrdr,maxrdr,tiny

tiny=1.0e-5

!scale for plot (x,y)
bb(2)=1
bb(1)=real(nwv)
minrdr=minval(rdr)+tiny !make sure axis are different
maxrdr=maxval(rdr)-tiny 
bb(3)=real(minrdr-0.1*(maxrdr-minrdr)) !reverse axes
bb(4)=real(maxrdr+0.1*(maxrdr-minrdr))

call pgvport(0.15,0.95,0.20,0.95) !make room around the edges for labels (-1.8)
call pgsci(1) !plotting colour
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS",0.0,0,"BCNTS",0.0,0) !draw axes and tics
call pglabel("Bandpass Channel","r/R*","")

allocate(xp(nwv),yp(nwv))
do i=1,nwv
   xp(i)=real(i)
   yp(i)=real(rdr(i))
enddo

!call pgpt(nwv,xp,yp,-1)
call pgline(nwv,xp,yp)


return
end subroutine plotrdr
