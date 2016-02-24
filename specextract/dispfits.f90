program dispfits
use precision
implicit none
integer :: iargc

if(iargc().lt.1)then
   write(0,*) "Usage: dispfits <Image1> [Image2] ... [ImageN]"
   write(0,*) "  <Image1> : First FITS image to display"
   write(0,*) "  [Image2] : Second FITS image to display (optional)"
   write(0,*) "  [ImageN] : Nth FITS image to display (optional)"
   stop
endif

end program dispfits
