# Makefile options compliant with our own installation of MUMPS 5.4.0 on irene

###########################################
# Environnement to load
###########################################
#source $ALL_CCCWORKDIR/fftw3.mumps.5.4.0.intel.20.rh8.module
# or
# module load feature/mkl/multi-threaded
# module load intel/20.0.0 mpi/openmpi/4.1.4 mkl/20.0.0
# module load blas/mkl/20.0.0 lapack/mkl/20.0.0 scalapack/mkl/20.0.0
# module load ptscotch/6.1.0 parmetis/4.0.3 metis/5.1.0
# module load fftw3/mkl/20.0.0
# module load subversion/1.14.0

RM = /bin/rm -f

CC =  mpicc 
FC =  mpif90 
FL =  mpif90 

F77 = mpif90 
F90 = mpif90 
ICC = mpicc 

AR = ar vr
RANLIB  = echo

#DOUBLE VERSUS SINGLE PRECISION LU FACTORIZATION:
#   -D_sp                   (single precision)
#   -D_dp                   (double precision)
EXT=
OPT_PRE = -cpp $(EXT)

#OPTIONS COMPILE
 
OPTF =  -O3 -assume byterecl -qopenmp $(OPT_PRE)
OPTC =  -O3 -assume byterecl -qopenmp $(OPT_PRE)
OPTL =  -O3 -assume byterecl -qopenmp $(OPT_PRE)

#DEBUG VERBEUX
#OPTF =  -assume byterecl -qopenmp -g -warn all $(OPT_PRE)
#OPTC =  -assume byterecl -qopenmp -g -warn all $(OPT_PRE)
#OPTL =  -assume byterecl -qopenmp -g -warn all $(OPT_PRE)

#DEBUG AGRESSIF
#OPTF =  -assume byterecl -qopenmp -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created $(OPT_PRE)
#OPTC =  -assume byterecl -qopenmp -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created $(OPT_PRE)
#OPTL =  -assume byterecl -qopenmp -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created $(OPT_PRE)

#DEBUG LIGHT
#OPTF =  -O2 -g -assume byterecl -qopenmp $(OPT_PRE)
#OPTC =  -O2 -g -assume byterecl -qopenmp $(OPT_PRE)
#OPTL =  -O2 -g -assume byterecl -qopenmp $(OPT_PRE)


###########################################
# OPTIMIZATION LIBRARY
###########################################
OPTIM=/ccc/work/cont003/gen10596/gen10596/svn/OPTIMIZATION_rh8
IOPTIM=$(OPTIM)/COMMON/include
LOPTIM=-L$(OPTIM)/lib -lSEISCOPE_OPTIM


########################################################################
# MUMPS LIBRARY
########################################################################
MUMPS_ROOT=/ccc/work/cont003/gen10596/gen10596/MUMPS_Builds/MUMPS_snapshot20221004_O2_g_intel20_rh8

#Begin orderings

SCOTCHDIR  = ${PTSCOTCH_ROOT}
ISCOTCH    = -I$(SCOTCHDIR)/include
LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr -lscotch  -lscotcherr -lscotchmetis

LPORDDIR = $(MUMPS_ROOT)/PORD/lib/
IPORD    = -I$(MUMPS_ROOT)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR = ${METIS_ROOT}/lib
IMETISDIR = -I${METIS_ROOT}/include
LPARMETISDIR = ${PARMETIS_ROOT}/lib
IPARMETISDIR = -I${PARMETIS_ROOT}/include

LMETIS = -L$(LPARMETISDIR) -lparmetis -L$(LMETISDIR) -lmetis
IMETIS = $(IMETISDIR) $(IPARMETISDIR)

#ORDERINGSF = -Dpord
ORDERINGSF  = -Dmetis -Dpord -Dparmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LPORD) $(LMETIS) $(LSCOTCH)
IORDERINGS = $(IPORD) $(IMETIS) $(ISCOTCH)

#End orderings

IMUMPS=$(MUMPS_ROOT)/include 
LMUMPS=-L$(MUMPS_ROOT)/lib -lsmumps -lmumps_common -lpord


########################################################################
# FFTW3 LIBRARY
########################################################################
FFTW3=${FFTW3_ROOT}
# or link with own libraries compiled from wind_codes/share/FFTW sources
#FFTW3=../../../../share/FFTW/fftw-3.3.8
IFFTW3=-I${FFTW3_INCDIR}
LFFTW3=-L${FFTW3_LIBDIR} -lfftw3


##############################################################################
# MKL LIBRARY
##############################################################################
#LMKL=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5  -lpthread
MKLDIR = $(MKLROOT)/lib/intel64
LMKL = -L$(MKLDIR) -lmkl_scalapack_lp64  -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -Wl,--end-group -qopenmp -lpthread


########################################################################
INCPAR=-I${IFFTW3} -I$(IOPTIM) -I$(IMUMPS) $(IORDERINGS) -I../inc
LIBPAR=${LFFTW3} $(LOPTIM) $(LMUMPS) $(LORDERINGS) $(LMKL)
########################################################################


INC = $(INCPAR)
LIB = $(LIBPAR)

