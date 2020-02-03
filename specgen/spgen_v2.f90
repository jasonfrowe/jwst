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
!random number vars
integer, dimension(3) :: now
integer :: seed
real(double) :: ran2,dumr
!local vars
integer :: i,noversample,nunit,filestatus,nmodeltype,iargc
real(double) :: xout, yout,rv,b
character(80) :: cline,modelfile

interface
	subroutine writefitsphdu(fileout,funit)
		use precision
     	implicit none
     	integer :: funit
     	character(200), dimension(3) :: fileout
	end subroutine writefitsphdu
end interface

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Command line arguments
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
if(iargc().lt.3)then
   write(0,*) "Usage: spgen <specmodel> <noversample> <planetmodel>"
   write(0,*) "   <specmodel> - Atlas-9 stellar model"
   write(0,*) " <noversample> - is new sampling for Kernel (must be > 0)"
   write(0,*) " <planetmodel> - name of planet model (A, rprs)"
   write(0,*) "           [b] - impact parameter - optional (must be > 0)"
   stop
endif

if(iargc().ge.4)then
   call getarg(4,cline)
   read(cline,*) b
   if(b.lt.0.0d0)then
      write(0,*) "b must be positive"
      stop
   endif
else
   !default impact parameter
   b=2.0!0.3589
endif

noversample=1 !now a commandline-parameter
!get oversampling from commandline
call getarg(2,cline)
read(cline,*) noversample !read in noversample
if(noversample.le.0)then
   write(0,*) "noversample must be greater than zero"
   stop
endif

!read in a model spectrum
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!read in a model spectrum
call getarg(1,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  Model Parameters 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!image dimensions
xout=2048  !dimensions for output image.
yout=256

!parameters that control the simulation
rv=0.0 !radial velocity shift (m/s)

!parameter controling modeltype
nmodeltype=2 !1=BT-Settl, 2=Atlas-9+NL limbdarkening

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!file naming
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!create FITS files and insert primary HDU for each output data product. 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 1 - simulation with no convolution,  Native resolution
! 2 - simulation with convolution, native resolution
! 3 - simulation with convolution, over-sampled resolution
call writefitsphdu(fileout(1),funit(1))
call writefitsphdu(fileout(2),funit(2))
call writefitsphdu(fileout(3),funit(3))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Random Number Initialization
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!Initialization of random number
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)




!close the FITS file
call closefits(funit(1))
call closefits(funit(2))
call closefits(funit(3))

end program spgenV2