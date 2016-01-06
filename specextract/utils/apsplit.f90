subroutine apsplit(naxes,Image,nlines,nTrace,solpsf,apfluxl,apfluxu)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer nlines,nTrace,i,j,k,n1,n2,nplot,nstart,nend,nlinesex
integer, dimension(2) :: naxes
real, allocatable, dimension(:) :: px,py
real(double) :: pos,width,apfac,avgsplitlength
real(double), dimension(:,:) :: Image,solpsf,apfluxl,apfluxu

apfac=2.0 !how big to make the flux aperture relative to Gaussian width

do k=1,nTrace
   !lower side of spectrum
   do i=1,nlines
      pos=solpsf(i,3+9*(k-1))+solpsf(i,9+9*(k-1))
      width=solpsf(i,4+9*(k-1))
      n1=max(1,int(pos-width*apfac-0.5))
      n2=min(naxes(2),int(pos+width*apfac+0.5))
      apfluxl(i,k)=Sum(Image(i,n1:n2))
!      write(0,*) "apfluxl:",apfluxl(i,k),n1,n2
   enddo

   !upper side of spectrum
   do i=1,nlines
      pos=solpsf(i,6+9*(k-1))+solpsf(i,9+9*(k-1))
      width=solpsf(i,4+9*(k-1))
      n1=max(1,int(pos-width*apfac-0.5))
      n2=min(naxes(2),int(pos+width*apfac+0.5))
      apfluxu(i,k)=Sum(Image(i,n1:n2))
!      write(0,*) "apfluxu:",apfluxu(i,k),n1,n2
   enddo

enddo


end subroutine apsplit
