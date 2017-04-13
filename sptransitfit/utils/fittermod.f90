module fittermod
   use precision, only: double
   implicit none
   integer, pointer :: nwv2,nobs2,nplanet2,npars2
   integer, dimension(:), pointer :: ntt2
   integer, dimension(:,:), pointer :: solrange2
   real(double), dimension(:), pointer :: sol2
   real(double), dimension(:,:), pointer ::  time2,flux2,ferr2,exptime2,&
    tobs2,omc2,solerr2

end module fittermod
