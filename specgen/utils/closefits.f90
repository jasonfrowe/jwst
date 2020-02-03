subroutine closefits(funit)
!simple routine that closes the active FITS file.
use precision
implicit none
integer :: funit,status

status=0

!close fits file
call ftclos(funit,status)
write(0,*) funit,status
call ftfiou(funit,status)
write(0,*) funit,status

return
end subroutine
