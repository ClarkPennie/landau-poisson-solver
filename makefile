# Directories
DIR=$(PWD)/
EXECDIR=$(DIR)
OBJDIR=$(DIR)
SRCDIR=$(DIR)

# GNU C compiler
MPICC=mpicxx 
CC=icc
CPP=icpc

#Intel MPI compiler for C++
MPICC=mpicxx

# Compiler flags: crashed when compiling with -O0
CFLAGS = -O2 -qopenmp  -I$(TACC_FFTW3_INC) -I$(TACC_MKL_INC)
FFTFLAGS = -L$(TACC_FFTW3_LIB) -lfftw3_threads -lfftw3 -lpthread -lm 
MKLFLAGS = -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread

# Command definition
RM=rm -f

sources_FPL = FPL_main.cpp 
objects_FPL= $(sources_FPL:.c=.o)

sources_LP = LP_ompi.cpp 
objects_LP= $(sources_LP:.c=.o)

sources_wt = WeightGenerator_mpi.cpp 
objects_wt= $(sources_wt:.c=.o)


FPL: $(objects_FPL)
	@echo "Building FPL solver"
	$(CPP) $(objects_FPL) $(CFLAGS)  -o fpl.out $(MKLFLAGS) $(FFTFLAGS)

LP: $(objects_LP)
	@echo "Building Landau-Poisson solver"
	@$(MPICC) $(CFLAGS) $(objects_LP) -o $(EXECDIR)LPsolver_nu005_TestNewFiles.out  $(FFTFLAGS) $(MKLFLAGS)

wts: $(objects_wt)
	@echo "Building mpi weights"
	@$(MPICC) -O2 -openmp $(objects_wt) -o $(EXECDIR)weight.out  

clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(EXECDIR)*.out

# icpc -O2 -openmp LP_main.cpp -I$TACC_FFTW3_INC -I$TACC_MKL_INC -L$TACC_FFTW3_LIB -lfftw3_threads -lfftw3 -lpthread -lm -Wl,-rpath,$TACC_MKL_LIB -L$TACC_MKL_LIB -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread