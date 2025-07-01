include ./Makefile.inc

EXEDIR = ../bin

.SUFFIXES: .f .o
.f.o:
	$(F77) $(OPTF) $(INC) -c $*.f

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(OPTF) $(INC) -c $*.f90

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -c $*.c

all: compile

compile: twri2d_mpi_v1

clean:
	$(RM) *.o *.mod *~
	
SUB_v1 = subcsrcal.f90 libfourier.o submodule.o kwsinc.o minresModule.o subget_parameter_gradient.o subproximity.o subadmm_algo.o \
         twri2d_mpi_v1.o subacqui.o subperturbation.o  subreadpar.o submodeling2.o subgradient.o readacqui.o \
         subresiduals.o submisfit.o subcor.o submodels.o sublinearminres.o \
         scalegradient.o subname.o para_range.o subresample_hicks.o subsmooth.o subfilmartin2d.o

twri2d_mpi_v1 : $(SUB_v1)
	$(F90) $(OPTF) $(SUB_v1) $(LIB) -o $(EXEDIR)/twri2d_mpi_v1
