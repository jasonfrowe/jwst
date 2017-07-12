# Makefile Template
F90 = ifort
LFLAGS =  -L/usr/lib -L/usr/local/lib -L/usr/X11/lib
XFLAGS = -lX11
PFLAGS = -lpng -lpgplot
CFLAGS = -lcfitsio
FFTFLAGS = -lfftw3
FFLAGS = -O3
BIN = ./bin/
UTILS = utils/

all: tilt apextract specextract calcldco spgen genoskernel

specextractsrc = specextract/

tiltincl = $(specextractsrc)precision.o $(specextractsrc)getfits.o $(specextractsrc)displayfits.o $(specextractsrc)rqsort.o $(specextractsrc)stdev.o $(specextractsrc)heatlut.o $(specextractsrc)apsplit.o $(specextractsrc)xcorr.o $(specextractsrc)spline.o $(specextractsrc)lfit.o $(specextractsrc)fitline.o $(specextractsrc)polyfilter.o
tilt: $(specextractsrc)apextract.f90 $(tiltincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ $(specextractsrc)tilt.f90 $(tiltincl) $(LFLAGS) $(XFLAGS) $(PFLAGS) $(CFLAGS) $(FFTFLAGS)

apextractincl = $(specextractsrc)precision.o $(specextractsrc)getfits.o $(specextractsrc)displayfits.o $(specextractsrc)rqsort.o $(specextractsrc)stdev.o $(specextractsrc)heatlut.o $(specextractsrc)apsplit.o $(specextractsrc)xcorr.o $(specextractsrc)spline.o $(specextractsrc)lfit.o $(specextractsrc)fitline.o $(specextractsrc)polyfilter.o
apextract: $(specextractsrc)apextract.f90 $(apextractincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ $(specextractsrc)apextract.f90 $(apextractincl) $(LFLAGS) $(XFLAGS) $(PFLAGS) $(CFLAGS) $(FFTFLAGS)

specextractincl = $(specextractsrc)precision.o $(specextractsrc)fittingmod.o $(specextractsrc)getfits.o $(specextractsrc)displayfits.o $(specextractsrc)rqsort.o $(specextractsrc)stdev.o $(specextractsrc)trace.o $(specextractsrc)apflux.o $(specextractsrc)heatlut.o $(specextractsrc)loadpsfinit.o $(specextractsrc)triplegaussian.o $(specextractsrc)psfmodel1d.o $(specextractsrc)modelline.o $(specextractsrc)minpack.o $(specextractsrc)getsky.o
specextract: $(specextractsrc)specextract.f90 $(specextractincl)
	$(F90) $(FFLAGS) -o $(BIN)$@ $(specextractsrc)specextract.f90 $(specextractincl) $(LFLAGS) $(XFLAGS) $(PFLAGS) $(CFLAGS)
	
calcldcosrc = limbdark/
calcldcoincl = $(calcldcosrc)precision.o $(calcldcosrc)minpack.o $(calcldcosrc)fittingmod.o
calcldco: $(calcldcosrc)calcldco.f90 $(calcldcoincl)
	$(F90) $(FFLAGS) -o $(BIN)$@ $(calcldcosrc)calcldco.f90 $(calcldcoincl) $(LFLAGS) $(XFLAGS) $(PFLAGS)

spgensrc = specgen/
spgenincl = $(spgensrc)precision.o $(spgensrc)response.o $(spgensrc)readmodel.o $(spgensrc)writefits.o $(spgensrc)deletefile.o $(spgensrc)binmodel.o $(spgensrc)rqsort.o $(spgensrc)addflux2pix.o $(spgensrc)readkernels.o $(spgensrc)getfits.o $(spgensrc)genkernel.o $(spgensrc)spline.o $(spgensrc)convolveft.o $(spgensrc)ovrwrt.o $(spgensrc)readresponse.o $(spgensrc)readheader.o
spgen: $(spgensrc)spgen.f90 $(spgenincl)
	$(F90) $(FFLAGS) -o $(BIN)$@ $(spgensrc)spgen.f90 $(spgenincl) $(LFLAGS) $(CFLAGS) $(FFTFLAGS)

genoskernelincl = $(spgensrc)precision.o $(spgensrc)readkernels.o $(spgensrc)getfits.o $(spgensrc)spline.o $(spgensrc)writefits.o $(spgensrc)deletefile.o
genoskernel: $(spgensrc)genoskernel.f90 $(genoskernelincl)
	$(F90) $(FFLAGS) -o $(BIN)$@ $(spgensrc)genoskernel.f90 $(genoskernelincl) $(LFLAGS) $(CFLAGS)

	
#%.o: $(src)$(UTILS)%.f90
#	$(F90) -c $(FFLAGS) $< -o $@
#%.o: $(src)$(UTILS)%.f
#	$(F90) -c $(FFLAGS) $< -o $@

#Libs for SPECEXTRACT

$(specextractsrc)polyfilter.o: $(specextractsrc)$(UTILS)polyfilter.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)polyfilter.f90
$(specextractsrc)fitline.o: $(specextractsrc)$(UTILS)fitline.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)fitline.f90
$(specextractsrc)lfit.o: $(specextractsrc)$(UTILS)lfit.f
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)lfit.f
$(specextractsrc)spline.o: $(specextractsrc)$(UTILS)spline.f
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)spline.f
$(specextractsrc)xcorr.o: $(specextractsrc)$(UTILS)xcorr.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)xcorr.f90
$(specextractsrc)apsplit.o: $(specextractsrc)$(UTILS)apsplit.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)apsplit.f90
$(specextractsrc)getsky.o: $(specextractsrc)$(UTILS)getsky.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)getsky.f90
$(specextractsrc)fittingmod.o: $(specextractsrc)$(UTILS)fittingmod.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)fittingmod.f90
$(specextractsrc)loadpsfinit.o: $(specextractsrc)$(UTILS)loadpsfinit.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)loadpsfinit.f90
$(specextractsrc)triplegaussian.o: $(specextractsrc)$(UTILS)triplegaussian.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)triplegaussian.f90
$(specextractsrc)psfmodel1d.o: $(specextractsrc)$(UTILS)psfmodel1d.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)psfmodel1d.f90
$(specextractsrc)modelline.o: $(specextractsrc)$(UTILS)modelline.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)modelline.f90
$(specextractsrc)minpack.o: $(specextractsrc)$(UTILS)minpack.f
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)minpack.f
$(specextractsrc)heatlut.o: $(specextractsrc)$(UTILS)heatlut.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)heatlut.f90
$(specextractsrc)apflux.o: $(specextractsrc)$(UTILS)apflux.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)apflux.f90 
$(specextractsrc)trace.o: $(specextractsrc)$(UTILS)trace.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)trace.f90 
$(specextractsrc)displayfits.o: $(specextractsrc)$(UTILS)displayfits.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)displayfits.f90 
$(specextractsrc)stdev.o: $(specextractsrc)$(UTILS)stdev.f
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)stdev.f
$(specextractsrc)precision.o: $(specextractsrc)$(UTILS)precision.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)precision.f90
$(specextractsrc)rqsort.o: $(specextractsrc)$(UTILS)rqsort.f
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)rqsort.f
$(specextractsrc)getfits.o: $(specextractsrc)$(UTILS)getfits.f90
	cd $(specextractsrc) ; $(F90) -c $(FFLAGS) $(UTILS)getfits.f90

$(calcldcosrc)precision.o: $(calcldcosrc)$(UTILS)precision.f90
	cd $(calcldcosrc) ; $(F90) -c $(FFLAGS) $(UTILS)precision.f90
$(calcldcosrc)minpack.o: $(calcldcosrc)$(UTILS)minpack.f
	cd $(calcldcosrc) ; $(F90) -c $(FFLAGS) $(UTILS)minpack.f
$(calcldcosrc)fittingmod.o: $(calcldcosrc)$(UTILS)fittingmod.f90
	cd $(calcldcosrc) ; $(F90) -c $(FFLAGS) $(UTILS)fittingmod.f90
	
$(spgensrc)response.o: $(spgensrc)$(UTILS)response.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)response.f90
$(spgensrc)readheader.o: $(spgensrc)$(UTILS)readheader.f
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)readheader.f
$(spgensrc)readresponse.o: $(spgensrc)$(UTILS)readresponse.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)readresponse.f90
$(spgensrc)precision.o: $(spgensrc)$(UTILS)precision.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)precision.f90
$(spgensrc)readmodel.o: $(spgensrc)$(UTILS)readmodel.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)readmodel.f90
$(spgensrc)writefits.o: $(spgensrc)$(UTILS)writefits.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)writefits.f90
$(spgensrc)deletefile.o: $(spgensrc)$(UTILS)deletefile.f
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)deletefile.f
$(spgensrc)binmodel.o: $(spgensrc)$(UTILS)binmodel.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)binmodel.f90
$(spgensrc)rqsort.o: $(spgensrc)$(UTILS)rqsort.f
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)rqsort.f
$(spgensrc)addflux2pix.o: $(spgensrc)$(UTILS)addflux2pix.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)addflux2pix.f90
$(spgensrc)readkernels.o: $(spgensrc)$(UTILS)readkernels.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)readkernels.f90
$(spgensrc)getfits.o: $(spgensrc)$(UTILS)getfits.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)getfits.f90
$(spgensrc)genkernel.o: $(spgensrc)$(UTILS)genkernel.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)genkernel.f90
$(spgensrc)spline.o: $(spgensrc)$(UTILS)spline.f
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)spline.f
$(spgensrc)convolve.o: $(spgensrc)$(UTILS)convolve.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)convolve.f90
$(spgensrc)convolveft.o: $(spgensrc)$(UTILS)convolveft.f90
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)convolveft.f90
$(spgensrc)ovrwrt.o: $(spgensrc)$(UTILS)ovrwrt.f
	cd $(spgensrc) ; $(F90) -c $(FFLAGS) $(UTILS)ovrwrt.f
	
# Removing object files
.PHONY : clean
clean :
	rm $(specextractsrc)*.o
	rm $(specextractsrc)*.mod
	rm $(calcldcosrc)*.o
	rm $(calcldcosrc)*.mod
	rm $(spgensrc)*.o
	rm $(spgensrc)*.mod
	
