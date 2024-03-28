# Configuration file for jorek
# model directory
MODEL = model502

# Fortran compiler 
FC     = mpif90
CC     = mpicc 
CXX    = mpicxx 

FFLAGS_OMP       = -fopenmp 
FFLAGS           = -fdefault-real-8 -fdefault-double-8 -DUSE_BLOCK -DFUNNELED -DWORLDWAR2 #-DCOMPARE_ELEMENT_MATRIX
FFLAGS          := $(FFLAGS) -O3 -msse2 -march=native

# --- DEBUGGING
#   --- Debug symbols for debuggers
DEBUGFLAGS  = #-g -p
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -traceback

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

#FFLAGS_NO_OMP    = $(FFLAGS) $(DEBUGFLAGS)


# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_PASTIX6 = 0
USE_MUMPS = 1
USE_WSMP   = 0
USE_FFTW = 1
USE_HDF5 = 1
USE_BOOST = 1
USE_STD_BESSELK = 0
USE_TASKLOOP = 1

#Scotch library
SCOTCH_HOME  = $(SCOTCH_ROOT)
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptesmumps -lptscotch -lesmumps -lscotch -lptscotcherrexit -lptscotcherr -lscotcherrexit -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# METIS
METIS_HOME = $(METIS_ROOT)
LIB_METIS = -L$(METIS_HOME)/lib -lmetis
INC_METIS = -I$(METIS_HOME)/include

# PARMETIS
PARMETIS_HOME = $(PARMETIS_ROOT)
LIB_PARMETIS = -L$(PARMETIS_HOME)/lib -lparmetis
INC_PARMETIS = -I$(PARMETIS_HOME)/include

# MUMPS
MUMPS_HOME = $(MUMPS_ROOT)
LIB_MUMPS  = -L$(MUMPS_HOME)/lib -lzmumps -ldmumps -lmumps_common $(LIB_SCOTCH) $(LIB_PARMETIS) $(LIB_METIS) 
INC_MUMPS  = $(MUMPS_HOME)/include $(INC_SCOTCH) $(INC_PARMETIS) $(INC_METIS)
ORDLIB     = -L$(MUMPS_HOME)/PORD/lib -lpord

# PASTIX
LIB_PASTIX        = `$(PASTIX_ROOT)/pastix-conf --libs` -L$(MKLLIB) -lpciaccess
LIB_PASTIX_MURGE  = `$(PASTIX_ROOT)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS   = `$(PASTIX_ROOT)/pastix-conf --blas`
INC_PASTIX        = `$(PASTIX_ROOT)/pastix-conf --incs`

PASTIX6_HOME      = -I$(PASTIX_ROOT)
LIB_PASTIX6       = -L$(PASTIX6_HOME)/lib -lpastixf -lspmf -lpastix -lspm -lpastix_kernels $(LIB_SCOTCH)
LIB_PASTIX6_BLAS  = 
INC_PASTIX6       = $(PASTIX6_HOME)/include $(INC_SCOTCH)

# BOOST
BOOST_HOME = $(BOOST_ROOT)
LIB_BOOST  = -L$(BOOST_HOME)/lib -lboost_math_tr1
INC_BOOST  = -I$(BOOST_HOME)/include 

#HIPS
#LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
#INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#WSMP
#LIB_WSMP = -L$(WSMP_HOME) -lpwsmp64

# LIBFFTW
FFTW_HOME = $(FFTW_ROOT)
LIBFFTW   = -L$(FFTW_HOME)/lib -lfftw3
INC_FFTW  = -I$(FFTW_HOME)/include

# HDF5
HDF5INCLUDE = $(HDF5_ROOT)/include
HDF5LIB     = $(HDF5_ROOT)/lib/libhdf5hl_fortran.a \
	      $(HDF5_ROOT)/lib/libhdf5_hl.a \
	      $(HDF5_ROOT)/lib/libhdf5_fortran.a \
	      $(HDF5_ROOT)/lib/libhdf5.a -lz -ldl

LIBLAPACK   = -L$(MKLLIB)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_blas95_lp64 -lmkl_tbb_thread -lmkl_gnu_thread -lmkl_sycl -lmkl_gf_lp64 -lmkl_cdft_core -lmkl_core -lmkl_sequential -lmkl_rt
