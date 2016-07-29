import numpy as np 
import nirisstools as nt

def tracespec(scidatain,nline=800,ntrace=1):
    "Finds trace of spectral orders"
    
    isize=scidatain.shape  #need to check size of input image.  transpose if necessary.
    if(isize[1]>isize[0]):
        scidata=np.transpose(scidatain)
    else :
        scidata=np.copy(scidatain)
        
    naxes=np.zeros(2,dtype="int32") #store size of image.
    naxes[0]=int(scidata.shape[0])
    naxes[1]=int(scidata.shape[1])
    
    bpix=np.float(1.0e30) #identifying bad pixels
    
    dtrace=np.zeros(shape=(naxes[0],ntrace),order='F')
    bf=np.zeros(shape=(naxes[0],ntrace),order='F')
    nfit=int(1+9*ntrace)
    solpsf=np.zeros(shape=(naxes[0],nfit),order='F')
    posguess=np.ones(ntrace)*(-1)
    
    #print(posguess)
    nt.trace(naxes,scidata,bpix,nline,ntrace,dtrace,bf,solpsf,posguess)
    return dtrace, solpsf ;
    
def apertureflux(scidatain,dtrace,napin,nskyin):
    "Extract Flux using simple Aperture"
    
    isize=scidatain.shape  #need to check size of input image.  transpose if necessary.
    if(isize[1]>isize[0]):
        scidata=np.transpose(scidatain)
    else :
        scidata=np.copy(scidatain)
        
    naxes=np.zeros(2,dtype="int32") #store size of image.
    naxes[0]=int(scidata.shape[0])
    naxes[1]=int(scidata.shape[1])
    
    bpix=np.float(1.0e30) #identifying bad pixels
    
    ntrace=int(dtrace.shape[1])
    
    nap=int(napin)   #make sure we have integer values
    nsky=int(nskyin)
    
    flux=np.zeros(dtrace.shape[0])
    nt.apflux(naxes,scidata,bpix,ntrace,dtrace,nap,nsky,flux)
    
    return flux;
    
