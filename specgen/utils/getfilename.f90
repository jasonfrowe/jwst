subroutine getfilename(pid,onum,vnum,gnum,spseq,anumb,enum,detectorname,prodtype,fileout)
use precision
implicit none
!import vars
integer :: pid,onum,vnum,gnum,spseq,anumb,enum
character(8) :: detectorname,prodtype
character(200), dimension(3) :: fileout !name of FITS files
!local vars
character(200) :: basename

!create basename
write(basename,'(A2,I0.5,I0.3,I0.3,A1,I0.2,I1,I0.2,I0.5,A1,A,A1,A)') 'jw', &
	pid,onum,vnum,'_',gnum,spseq,anumb,enum,'_',detectorname,'_',prodtype

!There are three FITS files made. 
!1) unconvolved Image
!2) convolved Image
!3) convolved Image with oversampled resolution.

fileout(1) = trim(basename)//'.fits'
fileout(2) = trim(basename)//'_c.fits'
fileout(3) = trim(basename)//'_ovc.fits'

return
end subroutine getfilename
