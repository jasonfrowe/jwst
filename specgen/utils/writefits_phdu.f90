subroutine writefitsphdu(fileout,funit)
!Open, create FITS file and insert PHDU with primary header.
!The primary HDU does not contain image data. 
use precision
implicit none
!import vars
integer :: funit !Unit used for FITS writing (In/Out)
character(200), dimension(3) :: fileout
!local vars
integer :: status,blocksize,bitpix,naxis
real(double), dimension(:), allocatable :: naxes
logical :: simple,extend


status=0
!if file already exists.. delete it.
call deletefile(fileout,status)
!get a unit number
call ftgiou(funit,status)

!Create the new empty FITS file.  The blocksize parameter is a
!historical artifact and the value is ignored by FITSIO.
blocksize=1
status=0
call ftinit(funit,fileout,blocksize,status)
if(status.ne.0)then
   write(0,*) "Status: ",status
   write(0,*) "Critial Error open FITS for writing"
   write(0,'(A80)') fileout
endif

!Initialize parameters about the FITS image.
!BITPIX = 16 means that the image pixels will consist of 16-bit
!integers.  The size of the image is given by the NAXES values.
!The EXTEND = TRUE parameter indicates that the FITS file
!may contain extensions following the primary array.
simple=.true.
bitpix=-32
naxis=1
allocate(naxes(1))
naxes(1)=1
extend=.true.

!Write the required header keywords to the file
call ftphpr(funit,simple,bitpix,naxis,naxes,0,1,extend,status)

!write(6,*) "ftprec:",status

!Adding required records for FITS file
write(6,*) 'JWST cards..'
call ftpkys(funit,'DATE-OBS','01/22/2020','/ [DD/MM/YYYY] Date of observation',status)
call ftpkyj(funit,'NRSTSTRT',1,'/ the number of resets at the start of the exposure',status)
call ftpkys(funit,'NRESETS',1,'/ the number of resets between integrations',status)
call ftpkys(funit,'DATE','2019-12-05T11:09:28.097', '/ [yyyy-mm-ddThh:mm:ss.ss] UTC date file cre',status)
!FILENAME= 'smalloutput.fits'   / Name of the file                               
call ftpkys(funit,'DATAMODL','RampModel','/ Type of data model',status)                             
call ftpkys(funit,'TELESCOP','JWST    ','/ Telescope used to acquire the data',status)
!Observation identifiers                                            
call ftpkys(funit,'TIME-OBS','11:08:45','/ [hh:mm:ss.sss] UTC time at start of exposure',status)  
!Target information                              
call ftpkyd(funit,'TARG_RA',188.38685,5,'/ Target RA at mid time of exposure',status)
call ftpkyd(funit,'TARG_DEC',-10.14617305555556,15,'/ Target Dec at mid time of exposure',status)  
call ftpkys(funit,'SRCTYPE','POINT   ','/ Advised source type (point/extended)',status)           
!Instrument configuration information                                               
call ftpkys(funit,'INSTRUME','NIRISS  ','/ Instrument used to acquire the data',status)   
call ftpkys(funit,'DETECTOR','NIS     ','/ Name of detector used to acquire the data',status)
call ftpkys(funit,'FILTER','CLEAR   ','/ Name of the filter element used',status)       
call ftpkys(funit,'PUPIL','GR700XD ','/ Name of the pupil element used',status)             
!Exposure parameters                                          
call ftpkys(funit,'EXP_TYPE','NIS_SOSS','/ Type of data in the exposure',status)             
call ftpkys(funit,'READPATT','NISRAPID','/ Readout pattern',status)                                
call ftpkyj(funit,'NINTS',1,'/ Number of integrations in exposure',status)         
call ftpkyj(funit,'NGROUPS',1,'/ Number of groups in integration',status)         
call ftpkyj(funit,'NFRAMES',1,'/ Number of frames per group',status)            
call ftpkyj(funit,'GROUPGAP',0,'/ Number of frames dropped between groups',status)
call ftpkyd(funit,'TFRAME',5.491,3,'/ [s] Time between frames',status)                        
call ftpkyd(funit,'TGROUP',5.491,3,'/ [s] Time between groups',status)                  
call ftpkyd(funit,'DURATION',0.0,1,'/ [s] Total duration of exposure',status)                         
!Subarray parameters                                                                           
call ftpkys(funit,'SUBARRAY','SUBSTRIP256','/ Subarray used',status)                             
call ftpkyj(funit,'SUBSTRT1',1,'/ Starting pixel in axis 1 direction',status)            
call ftpkyj(funit,'SUBSTRT2',1793,'/ Starting pixel in axis 2 direction',status)             
call ftpkyj(funit,'SUBSIZE1',2048,'/ Number of pixels in axis 1 direction',status)           
call ftpkyj(funit,'SUBSIZE2',256,'/ Number of pixels in axis 2 direction',status)        
call ftpkyj(funit,'FASTAXIS',-2,'/ Fast readout axis direction',status)        
call ftpkyj(funit,'SLOWAXIS',-1,'/ Slow readout axis direction',status)


return
end subroutine writefitsphdu