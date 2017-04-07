#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
FPATH = /opt/intel/compilers_and_libraries_2017.2.163/mac/bin/intel64/
F90 = $(FPATH)ifort
F77 = $(FPATH)ifort
#compiling object file flags
OPT1 = -O3 -ipo
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
BIN = /Users/rowe/Documents/jwst/sptransitfit/
#utils source directory
UTILS = utils/

#Listing of programs to create.
all: transitfit8

transitfit8incl = precision.o readdata.o getfitpars.o sptransitmodel.o ttcor.o keplerian.o albedomod.o mandelagol.o occultquad.o fittransitmodel8.o estzpt.o getb.o rqsort.o
transitfit8: transitfit8.f90 $(transitfit8incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitfit8.f90 $(transitfit8incl)

#building object libraries
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
getb.o: $(UTILS)getb.f90
	$(F90) $(FFLAGS) $(UTILS)getb.f90
estzpt.o: $(UTILS)estzpt.f90
	$(F90) $(FFLAGS) $(UTILS)estzpt.f90
fittransitmodel8.o: $(UTILS)fittransitmodel8.f90
	$(F90) $(FFLAGS) $(UTILS)fittransitmodel8.f90
occultquad.o: $(UTILS)occultquad.f
	$(F90) $(FFLAGS) $(UTILS)occultquad.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F90) $(FFLAGS) $(UTILS)mandelagol.f
albedomod.o: $(UTILS)albedomod.f90
	$(F90) $(FFLAGS) $(UTILS)albedomod.f90
keplerian.o: $(UTILS)keplerian.f
	$(F90) $(FFLAGS) $(UTILS)keplerian.f
ttcor.o: $(UTILS)ttcor.f
	$(F90) $(FFLAGS) $(UTILS)ttcor.f
sptransitmodel.o: $(UTILS)sptransitmodel.f90
	$(F90) $(FFLAGS) $(UTILS)sptransitmodel.f90
getfitpars.o: $(UTILS)getfitpars.f90
	$(F90) $(FFLAGS) $(UTILS)getfitpars.f90
readdata.o: $(UTILS)readdata.f90
	$(F90) $(FFLAGS) $(UTILS)readdata.f90
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
