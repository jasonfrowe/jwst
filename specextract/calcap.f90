program calcap
!Estimated aperture to enclose 90,95 and 99% of flux.
use precision
implicit none
integer iargc

if(iargc().lt.1)then
   write(0,*) "Usage : calcap <tracepsf.dat>"
   write(0,*) "  <tracepsf.dat> : Output from spectextract"
   stop
endif

end program calcap
