subroutine plotres(nwv,nobs,res)
use precision
implicit none
!import vars
integer :: nwv,nobs,ncol
real(double),dimension(:,:) :: res
!local vars
integer :: i,ia1
real :: bb(4),maxyp,minyp,r,g,b
real, allocatable, dimension(:) :: xp,yp

allocate(xp(nobs),yp(nobs))

!scale for plot (x,y)
bb(1)=1
bb(2)=real(nobs)
minyp=real(minval(res))
maxyp=real(maxval(res))
bb(3)=minyp-0.1*(maxyp-minyp)
bb(4)=maxyp+0.1*(maxyp-minyp)

call pgvport(0.15,0.95,0.20,0.95) !make room around the edges for labels (-1.8)
call pgsci(1) !plotting colour
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS",0.0,0,"BCNTS",0.0,0) !draw axes and tics
call pglabel("Time","Residuals","")

!set up plot colours
ncol=64 !number of colours for display
call pgscr(15,0.0,0.3,0.2)
call pgsci(1)

!load in colourmap
do i=1,ncol
   call heatlut(i*4-3,r,g,b)
   call heatlut(i*4-2,r,g,b)
   call heatlut(i*4-1,r,g,b)
   call heatlut(i*4  ,r,g,b)
   CALL PGSCR(I+15, R, G, B)
enddo

do i=1,nobs
   xp(i)=real(i)
enddo

do i=1,nwv

   ia1= int(real(i-1)/real(nwv-1)*real(ncol-1))+16
   if(i.lt.1)   ia1=16
   if(i.gt.nwv) ia1=ncol+15
   call pgsci(ia1)

   yp=real(res(i,:))
   call pgpt(nobs,xp,yp,-1)

enddo

call pgsci(1) !return to default plotting colour

return
end subroutine plotres
