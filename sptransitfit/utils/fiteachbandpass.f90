subroutine fiteachbandpass(npars,nplanet,sol,solerr,solrange,nwv,       &
 nobs,time,flux,ferr,exptime,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: npars,nwv,nobs,nplanet
integer, dimension(:) :: ntt
integer, dimension(:,:) :: solrange
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,tobs,omc
!local vars
integer i,ii,k,inwv,npars_1,np,nwv_1
integer, allocatable, dimension(:,:) :: solrange_1
real(double), allocatable, dimension(:) :: sol_1
real(double), allocatable, dimension(:,:) :: time_1,flux_1,ferr_1,      &
 exptime_1,tn_1,solerr_1
character(80) :: line

!interfaces to allow assumed arrays
interface
   subroutine fittransitmodel8(npars,nplanet,sol,solerr,solrange,nwv, &
    nobs,time,flux,ferr,exptime,ntt,tobs,omc)
      use precision
      implicit none
      integer :: npars,nwv,nobs,nplanet
      integer, dimension(:) :: ntt
      integer, dimension(:,:) :: solrange
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: solerr,time,flux,ferr,exptime,    &
       tobs,omc
   end subroutine fittransitmodel8
end interface

!allocate array space for single bandpass photometry
allocate(time_1(1,nobs),flux_1(1,nobs),ferr_1(1,nobs),exptime_1(1,nobs))

!allocate arrays to hold solution for each single bandpass
npars_1=8+10*nplanet !number of parameters for whitelight model
allocate(sol_1(npars_1),solerr_1(npars_1,4),solrange_1(npars_1,2))

do inwv=1,nwv

   write(line,'(A19,I5)') "Fitting Bandpass # ",inwv
   call ovrwrt(line,2)

   !create weighted averages across all bandpasses
   time_1(1,:)=time(inwv,:) !time array for one bandpass
   flux_1(1,:)=flux(inwv,:) !flux array for one bandpass
   ferr_1(1,:)=ferr(inwv,:) !uncertainties for one bandpass
   exptime_1(1,:)=exptime(inwv,:) !exposure times for one bandpass

   do i=1,8 !loop over global parameters
      if(solrange(i,2)-solrange(i,1).gt.0)then
         k=inwv-1
      else
         k=0
      endif
      sol_1(i)=sol(solrange(i,1)+k)
      solerr_1(i,:)=solerr(solrange(i,1)+k,:)
      solrange_1(i,:)=i
   enddo
   !disable fit for rhostar
   solerr_1(1,1)=0.0d0

   !planet parameters
   do np=1,nplanet !loop over each planet in transit model
      do i=1,10    !loop over each parameter for a single planet
         ii=8+i+nplanet*(np-1)
         if(solrange(ii,2)-solrange(ii,1).gt.0)then
            k=inwv-1
         else
            k=0
         endif
         sol_1(ii)=sol(solrange(ii,1)+k)
         solerr_1(ii,:)=solerr(solrange(ii,1)+k,:)
         solrange_1(ii,:)=ii
         if((i.eq.1).or.(i.eq.2).or.(i.eq.3).or.(i.eq.5).or.(i.eq.6).or.&
          (i.eq.7))then
            solerr_1(ii,1)=0.0d0 !disable b,ecosw,esinw,K
         endif
      enddo
   enddo

   !do i=1,npars_1
   !   write(0,500) i,sol_1(i)
   !enddo

   nwv_1=1 !there is only one band pass for whitelight
   !Fit model to single bandpass
   call fittransitmodel8(npars_1,nplanet,sol_1,solerr_1,solrange_1,     &
      nwv_1,nobs,time_1,flux_1,ferr_1,exptime_1,ntt,tobs,omc)

   !copy solution for individual bandpass into master
   do i=1,npars_1
      if(solrange(i,2)-solrange(i,1).gt.0)then
         k=inwv-1
      else
         k=0
      endif
      sol(solrange(i,1)+k)=sol_1(i)
   enddo

   !do i=1,npars_1
   !   write(0,500) i,sol_1(i)
   !enddo
   !500 format(I2,1X,1PE17.10)
   !write(0,*) "Ready..."
   !read(5,*)

enddo
write(0,*) " "

return
end subroutine fiteachbandpass
