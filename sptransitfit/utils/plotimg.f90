subroutine plotimg(nwv,nobs,res)
use precision
implicit none
!import vars
integer :: nwv,nobs
real(double), dimension(:,:) :: res
!local vars
integer :: i,j,ncol,ia1
integer, allocatable, dimension(:,:) :: ia
real:: r,g,b,bb(4)
real(double) :: resmax,resmin
real(double), allocatable, dimension(:,:) :: resback

allocate(resback(nwv,nobs))
resback=res

!scale for plot (x,y)
bb(1)=1
bb(2)=real(nobs)
bb(3)=1
bb(4)=real(nwv)

call PGSCLP(0) !turn off clipping
call pgvport(0.15,0.95,0.20,1.95) !make room around the edges for labels
call pgsci(1)
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS",0.0,0,"BCNTS",0.0,0)

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

!res=abs(res)
!res=log10(res+1.0e-6)
!res=res*abs(resback)/res
!get max/min of plotting scale
resmax=maxval(res)
resmin=minval(res)


allocate(ia(nobs,nwv))
do i=1,nwv
   do j=1,nobs
      ia1=int((res(i,j)-resmin)/(resmax-resmin)*dble(NCOL-1))+16
      if(res(i,j).lt.resmin) ia1=16
      if(res(i,j).gt.resmax) ia1=ncol+15
      ia(j,i)=ia1
   enddo
enddo

call pgpixl(ia,nobs,nwv,1,nobs,1,nwv,bb(1),bb(2),bb(3),bb(4))
call pgbox("BCNTS",0.0,0,"BCNTS",0.0,0) !redraw boundary
call pglabel("Time","Bandpass","")

res=resback

return
end subroutine plotimg
