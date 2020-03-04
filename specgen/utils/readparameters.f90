subroutine readparameters(nunit,tstart,tend,exptime,deadtime,modelfile,nmodeltype,rvstar,vsini, &
   pmodelfile,nplanet,nplanetmax,sol,xout,yout,noversample,saturation,ngroup, &
   pid,onum,vnum,gnum,spseq,anumb,enum,enumos,detectorname,prodtype)
use precision
implicit none
!import vars
integer :: nmodeltype,nplanet,nplanetmax,xout,yout,noversample,ngroup,pid,onum, &
   vnum,gnum,spseq,anumb,enum,enumos,nunit
real(double) :: tstart,tend,exptime,deadtime,rvstar,vsini,saturation
real(double) :: sol(nplanetmax*8+1)
character(8) :: detectorname,prodtype
character(80) :: modelfile,pmodelfile(nplanetmax)
!local parameters
integer :: i,ic,filestatus,np
real(double) :: dvalue
character(200) :: command,keyword

!default parameters -- these will cause the program to end quickly
tstart=0.0d0 !start time (hours)
tend=0.0d0 !end time (hours)
exptime=0.0d0 !exposure time (s)
deadtime=0.0d0 !dead time (s)
sol(1)=1.0d0 !mean stellar density (cgs)
modelfile='null' !stellar spectrum file name
nmodeltype=2 !stellar spectrum type. 1=BT-Settl, 2=Atlas-9+NL limbdarkening
rvstar=0.0d0 !radial velocity of star (km/s)
vsini=0.0d0 !projected rotation of star (km/s)
pmodelfile='null' !file with Rp/Rs values 
nplanet=0 !number of planets -- default is no planets - you will get staronly sim.
sol(2)=0.0 !T0 (days)
sol(3)=1.0 !period (days)
sol(4)=0.0 !impact parameter 
sol(5)=0.0 !sqrt(e)sinw
sol(6)=0.0 !sqrt(e)cosw
sol(7)=0.0 !KRV (m/s)
sol(8)=0.0 !geometric albedo
sol(9)=0.0 !ellipsodial variations
xout=2048  !dispersion axis
yout=256   !spatial axis
noversample=1 !oversampling
saturation=65536.0d0 !saturation
ngroup=1 !samples up ramp
pid = 1 !programID
onum = 1 !observation number
vnum = 1 !visit number
gnum = 1 !group visit
spseq = 1 !parallel sequence. (1=prime, 2-5=parallel)
anumb = 1 !activity number
enum = 1 !exposure number
enumos = 1 !exposure number for oversampling
detectorname = 'NISRAPID' !convert this 
prodtype='cal'

!Default is having no transiting planets.
do i=1,nplanetmax
   pmodelfile(i)='null'
enddo
nplanet=0 !count number of planets.

i=1 !count line number
!Read in parameters..
do
   read(nunit,500,iostat=filestatus) command
   500 format(A200) !command much be contained within first 80 characters
   
   if((command(1:1).ne.'#').and.(command(1:1).ne.' '))then !ignore comment and blank lines
      ic=scan(command,' ')
      read(command(1:ic),*) keyword 
      write(0,*) 'Keyword: ',keyword(1:ic)
      select case(keyword)
      case('TSTART')
         read(command(ic+1:),*) dvalue
         tstart=dvalue
      case('TEND')
         read(command(ic+1:),*) dvalue
         tend=dvalue
      case('EXPTIME')
         read(command(ic+1:),*) dvalue
         exptime=dvalue
      case('DEADTIME')
         read(command(ic+1:),*) dvalue
         deadtime=dvalue
      case('RHOSTAR')
         read(command(ic+1:),*) dvalue
         sol(1)=dvalue
      case('STARMOD')
         read(command(ic+1:),*) modelfile
      case('STARTYPE')
         read(command(ic+1:),*) dvalue
         nmodeltype=int(dvalue)
         if(nmodeltype.ne.2)then
            write(0,*) "Error: Only ATLAS-9 models are current support (STARTYPE=2)"
            stop
         endif
      case('VSINI')
         read(command(ic+1:),*) dvalue
         vsini=dvalue
      case default
         !right now, this will pick up planet parameters
         write(0,*) 'Warning: Invalid KEYWORD: ', keyword(1:ic)
      end select

      !handle multiplanet systems.
      if(keyword(1:ic-1).eq.'RPRSFILE')then
         read(keyword(ic:ic),*) np !get planet number
         if((np.le.9).and.(np.gt.0))then
            nplanet=max(nplanet,np) !keep track of how many planets we have.
            read(keyword(1:ic-1),*) pmodelfile(np)
         else
            write(0,*) 'Error: Planet number is Invalid ',np 
         endif
      endif

   endif

   if(filestatus == 0) then
      i=i+1
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

!convert all times to days
tstart=tstart/24.0d0 !hours -> days
tend=tend/24.0d0     !hours -> days
exptime=exptime/86400.0   !seconds -> days
deadtime=deadtime/86400.0 !seconds -> days

!check if modelfile is null.  If so, exit, because there is nothing to do.

return
end

