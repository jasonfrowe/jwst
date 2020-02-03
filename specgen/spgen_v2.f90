program spgenV2
!Version 2 of Spec-generator.  This version does time-series and get FITS formating correct. 
!generates GR700-SOSS spectrum with 3-order traces + PSF + response
use precision
!use response
implicit none
integer, dimension(3) :: funit !number of FITS I/O
character(200), dimension(3) :: fileout !name of FITS files
!file name vars
integer :: pid,onum,vnum,gnum,spseq,anumb,enum
character(8) :: detectorname,prodtype
!local vars
integer :: i,iargc
real(double) :: xout, yout,rv

interface
	subroutine writefitsphdu(fileout,funit)
		use precision
     	implicit none
     	integer :: funit
     	character(200), dimension(3) :: fileout
	end subroutine writefitsphdu
end interface

if(iargc().lt.3)then
   write(0,*) "Usage: spgen <specmodel> <noversample> <planetmodel> [b] [time]"
   write(0,*) "   <specmodel> - Atlas-9 stellar model"
   write(0,*) " <noversample> - is new sampling for Kernel (must be > 0)"
   write(0,*) " <planetmodel> - name of planet model (A, rprs)"
   write(0,*) "           [b] - impact parameter - optional (must be > 0)"
   write(0,*) "        [time] - observation timestamp to embed in header"
   write(0,*) " [output.fits] - filename for output"
   stop
endif

!image dimensions
xout=2048  !dimensions for output image.
yout=256

!parameters that control the simulation
rv=0.0 !radial velocity shift (m/s)

!file naming
pid = 1 !programID
onum = 1 !observation number
vnum = 1 !visit number
gnum = 1 !group visit
spseq = 1 !parallel sequence. (1=prime, 2-5=parallel)
anumb = 1 !activity number
enum = 1 !exposure number
detectorname = 'NISRAPID' !convert this 
prodtype='cal'
call getfilename(pid,onum,vnum,gnum,spseq,anumb,enum,detectorname,prodtype,fileout)

!create FITS files and insert primary HDU for each data type. 
! 1 - simulation with no convolution,  Native resolution
! 2 - simulation with convolution, native resolution
! 3 - simulation with convolution, over-sampled resolution
call writefitsphdu(fileout(1),funit(1))
call writefitsphdu(fileout(2),funit(2))
call writefitsphdu(fileout(3),funit(3))





!close the FITS file
call closefits(funit(1))
call closefits(funit(2))
call closefits(funit(3))

end program spgenV2