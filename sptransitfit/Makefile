#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
#FPATH = /opt/intel/compilers_and_libraries_2017.2.163/mac/bin/intel64/
FPATH = 
F90 = $(FPATH)ifort
F77 = $(FPATH)ifort
#compiling object file flags
OPT1 = -O3 -mkl -qopenmp -heap-arrays 
#OPT1 = -O0 -g -CB -mkl
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
X11DIR = /usr/X11/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib
#LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/jwst/sptransitfit/
#utils source directory
UTILS = utils/

#Listing of programs to create.
all: transitfit8 transitplot8 transitmcmc8

transitmcmc8incl = precision.o readdata.o getfitpars.o ttcor.o sptransitmodel.o \
 keplerian.o albedomod.o mandelagol.o occultquad.o ran2.o initran.o mhgmcmc.o \
 loglikelihood.o
transitmcmc8: transitmcmc8.f90 $(transitmcmc8incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitmcmc8.f90 $(transitmcmc8incl)

transitplot8incl = precision.o readdata.o getfitpars.o ttcor.o sptransitmodel.o \
 keplerian.o albedomod.o mandelagol.o occultquad.o heatlut.o plotimg.o plotrdr.o \
 plotres.o
transitplot8: transitplot8.f90 $(transitplot8incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitplot8.f90 $(transitplot8incl) $(LIBS)

transitfit8incl = precision.o fittermod.o readdata.o getfitpars.o \
 sptransitmodel.o ttcor.o keplerian.o albedomod.o mandelagol.o occultquad.o \
 fittransitmodel8.o fittransitmodel8v2.o  estzpt.o getb.o rqsort.o lbfgsb.o \
 linpack.o timer.o loglikelihood.o gradient.o exportfitpars.o gaussj.o covsrt.o \
 heatlut.o plotimg.o fitwhitelighttransit.o fiteachbandpass.o ovrwrt.o
transitfit8: transitfit8.f90 $(transitfit8incl)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitfit8.f90 $(transitfit8incl) $(LIBS)

#building object libraries
%.o : $(UTILS)%.f90
	$(F90) $(FFLAGS) -o $@ $<

%.o : $(UTILS)%.f
	$(F90) $(FFLAGS) -o $@ $<

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
