#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo
#OPT1 = -O0 -g -CB
#OPT1 = -O3
OPT2 =
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
all: spgen

genoskernelincl = precision.o readkernels.o getfits.o spline.o writefits.o deletefile.o
genoskernel: genoskernel.f90 $(genoskernelincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ genoskernel.f90 $(genoskernelincl) $(LIBS2)

spgenincl = precision.o response.o readmodel.o writefits.o deletefile.o binmodel.o rqsort.o addflux2pix.o readkernels.o getfits.o genkernel.o spline.o convolve.o ovrwrt.o readresponse.o readheader.o 
spgen: spgen.f90 $(spgenincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ spgen.f90 $(spgenincl) $(LIBS2) 


#building object libraries
response.o: $(UTILS)response.f90
	$(F90) $(FFLAGS) $(UTILS)response.f90
readheader.o: $(UTILS)readheader.f
	$(F90) $(FFLAGS) $(UTILS)readheader.f
readresponse.o: $(UTILS)readresponse.f90
	$(F90) $(FFLAGS) $(UTILS)readresponse.f90
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
readmodel.o: $(UTILS)readmodel.f90
	$(F90) $(FFLAGS) $(UTILS)readmodel.f90
writefits.o: $(UTILS)writefits.f90
	$(F90) $(FFLAGS) $(UTILS)writefits.f90
deletefile.o: $(UTILS)deletefile.f
	$(F90) $(FFLAGS) $(UTILS)deletefile.f
binmodel.o: $(UTILS)binmodel.f90
	$(F90) $(FFLAGS) $(UTILS)binmodel.f90
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
addflux2pix.o: $(UTILS)addflux2pix.f90
	$(F90) $(FFLAGS) $(UTILS)addflux2pix.f90
readkernels.o: $(UTILS)readkernels.f90
	$(F90) $(FFLAGS) $(UTILS)readkernels.f90
getfits.o: $(UTILS)getfits.f90
	$(F90) $(FFLAGS) $(UTILS)getfits.f90
genkernel.o: $(UTILS)genkernel.f90
	$(F90) $(FFLAGS) $(UTILS)genkernel.f90
spline.o: $(UTILS)spline.f
	$(F90) $(FFLAGS) $(UTILS)spline.f
convolve.o: $(UTILS)convolve.f90
	$(F90) $(FFLAGS) $(UTILS)convolve.f90
ovrwrt.o: $(UTILS)ovrwrt.f
	$(F90) $(FFLAGS) $(UTILS)ovrwrt.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
