#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#compiling object file flags
OPT1 = -O3 -parallel -ipo
#OPT1 = -O0 -g -CB
#OPT1 = -O3 -parallel -ipo
OPT2 = -heap-arrays 0 
#OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB
#fitsio libraries
FITSIODIR = /usr/local/lib
#Pgplot plot libraries
PGPLOTDIR = /usr/local/lib
#X11 libraries
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib
LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/jwst/specgen/
#utils source directory
UTILS = utils/

#Listing of programs to create.
all: specextract apextract tilt imintg

calcapincl = precision.o
calcap: calcap.f90 $(calcapincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ calcap.f90 $(calcapincl)

immatchincl = precision.o getfits.o displayfits.o heatlut.o stdev.o rqsort.o writefits.o deletefile.o
immatch: immatch.f90 $(immatchincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ immatch.f90 $(immatchincl) $(LIBS2)

imintgincl = precision.o readcfits.o
imintg: imintg.f90 $(imintgincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ imintg.f90 $(imintgincl) -L$(FITSIODIR) -lcfitsio

apextractincl = precision.o getfits.o displayfits.o rqsort.o stdev.o heatlut.o apsplit.o xcorr.o spline.o lfit.o fitline.o polyfilter.o
apextract: apextract.f90 $(apextractincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ apextract.f90 $(apextractincl) $(LIBS2) -lfftw3

tiltincl = precision.o getfits.o displayfits.o rqsort.o stdev.o heatlut.o apsplit.o xcorr.o spline.o lfit.o fitline.o polyfilter.o
tilt: tilt.f90 $(apextractincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ tilt.f90 $(tiltincl) $(LIBS2) -lfftw3

specextractincl = precision.o fittingmod.o getfits.o displayfits.o rqsort.o stdev.o trace.o apflux.o heatlut.o loadpsfinit.o triplegaussian.o psfmodel1d.o modelline.o minpack.o getsky.o
specextract: specextract.f90 $(specextractincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ specextract.f90 $(specextractincl) $(LIBS2)

#building object libraries
readcfits.o: $(UTILS)readcfits.f90
	$(F90) $(FFLAGS) $(UTILS)readcfits.f90
polyfilter.o: $(UTILS)polyfilter.f90
	$(F90) $(FFLAGS) $(UTILS)polyfilter.f90
fitline.o: $(UTILS)fitline.f90
	$(F90) $(FFLAGS) $(UTILS)fitline.f90
lfit.o: $(UTILS)lfit.f
	$(F90) $(FFLAGS) $(UTILS)lfit.f
spline.o: $(UTILS)spline.f
	$(F90) $(FFLAGS) $(UTILS)spline.f
xcorr.o: $(UTILS)xcorr.f90
	$(F90) $(FFLAGS) $(UTILS)xcorr.f90
apsplit.o: $(UTILS)apsplit.f90
	$(F90) $(FFLAGS) $(UTILS)apsplit.f90
getsky.o: $(UTILS)getsky.f90
	$(F90) $(FFLAGS) $(UTILS)getsky.f90
fittingmod.o: $(UTILS)fittingmod.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmod.f90
modelline.o: $(UTILS)modelline.f90
	$(F90) $(FFLAGS) $(UTILS)modelline.f90
minpack.o: $(UTILS)minpack.f
	$(F90) $(FFLAGS) $(UTILS)minpack.f
psfmodel1d.o: $(UTILS)psfmodel1d.f90
	$(F90) $(FFLAGS) $(UTILS)psfmodel1d.f90
loadpsfinit.o: $(UTILS)loadpsfinit.f90
	$(F90) $(FFLAGS) $(UTILS)loadpsfinit.f90
triplegaussian.o: $(UTILS)triplegaussian.f90
	$(F90) $(FFLAGS) $(UTILS)triplegaussian.f90
heatlut.o: $(UTILS)heatlut.f90
	$(F90) $(FFLAGS) $(UTILS)heatlut.f90
apflux.o: $(UTILS)apflux.f90
	$(F90) $(FFLAGS) $(UTILS)apflux.f90
trace.o: $(UTILS)trace.f90
	$(F90) $(FFLAGS) $(UTILS)trace.f90
displayfits.o: $(UTILS)displayfits.f90
	$(F90) $(FFLAGS) $(UTILS)displayfits.f90
stdev.o: $(UTILS)stdev.f
	$(F90) $(FFLAGS) $(UTILS)stdev.f
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
#readmodel.o: $(UTILS)readmodel.f90
#	$(F90) $(FFLAGS) $(UTILS)readmodel.f90
writefits.o: $(UTILS)writefits.f90
	$(F90) $(FFLAGS) $(UTILS)writefits.f90
deletefile.o: $(UTILS)deletefile.f
	$(F90) $(FFLAGS) $(UTILS)deletefile.f
#binmodel.o: $(UTILS)binmodel.f90
#	$(F90) $(FFLAGS) $(UTILS)binmodel.f90
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
#addflux2pix.o: $(UTILS)addflux2pix.f90
#	$(F90) $(FFLAGS) $(UTILS)addflux2pix.f90
#readkernels.o: $(UTILS)readkernels.f90
#	$(F90) $(FFLAGS) $(UTILS)readkernels.f90
getfits.o: $(UTILS)getfits.f90
	$(F90) $(FFLAGS) $(UTILS)getfits.f90
#genkernel.o: $(UTILS)genkernel.f90
#	$(F90) $(FFLAGS) $(UTILS)genkernel.f90
#spline.o: $(UTILS)spline.f
#	$(F90) $(FFLAGS) $(UTILS)spline.f
#convolve.o: $(UTILS)convolve.f90
#	$(F90) $(FFLAGS) $(UTILS)convolve.f90
#ovrwrt.o: $(UTILS)ovrwrt.f
#	$(F90) $(FFLAGS) $(UTILS)ovrwrt.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
