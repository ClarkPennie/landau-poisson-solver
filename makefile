# Directories
DIR=$(PWD)
EXECDIR:=$(DIR)/bin
OBJDIR:=$(DIR)/build
SRCDIR:=$(DIR)/source

# Files
EXEC :=  LPsolver_nu005_TestNewDoping.out 
SRC  :=  $(wildcard $(SRCDIR)/*.cpp) 
OBJ  :=  $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRC))

# Intel C compiler
CC=icc

# Intel C++ compiler
CPP=icpc

# Intel MPI compiler for C++
MPICC=mpicxx -g -p

# Compiler flags: crashed when compiling with -O0
CFLAGS = -O2 -qopenmp  -I$(TACC_FFTW3_INC) -I$(TACC_MKL_INC) -I$(SRCDIR)
FFTFLAGS = -L$(TACC_FFTW3_LIB) -lfftw3_threads -lfftw3 -lpthread -lm 
FFTINC = -I$(TACC_FFTW3_INC)
MKLFLAGS = -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread

# Command definition
RM=rm -f

sources_FPL = FPL_main.cpp 
objects_FPL= $(sources_FPL:.c=.o)

sources_LP = $(SRCDIR)/LP_ompi.cpp 
objects_LP= $(sources_LP:.c=.o)

sources_wt = WeightGenerator_mpi.cpp 
objects_wt= $(sources_wt:.c=.o)


FPL: $(objects_FPL)
	@echo "Building FPL solver"
	$(CPP) $(objects_FPL) $(CFLAGS)  -o fpl.out $(MKLFLAGS) $(FFTFLAGS)

LP: $(OBJ)
	@echo "Building Landau-Poisson solver"
	$(MPICC) $(CFLAGS) -o $(EXECDIR)/$(EXEC) $^	$(FFTFLAGS) $(MKLFLAGS)
	
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(MPICC) $(CFLAGS) $(FFTINC) -c -o $@ $<
	
$(OBJDIR)/LP_ompi.o: $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h $(SRCDIR)/collisionRoutines_1.h $(SRCDIR)/conservationRoutines.h $(SRCDIR)/EntropyCalculations.h $(SRCDIR)/EquilibriumSolution.h $(SRCDIR)/MarginalCreation.h $(SRCDIR)/MomentCalculations.h $(SRCDIR)/NegativityChecks.h $(SRCDIR)/FieldCalculations.h $(SRCDIR)/SetInit_1.h 
$(OBJDIR)/advection_1.o: $(SRCDIR)/advection_1.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/FieldCalculations.h
$(OBJDIR)/collisionRoutines_1.o: $(SRCDIR)/collisionRoutines_1.h $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h #$(SRCDIR)/ThreadPriv.h
$(OBJDIR)/conservationRoutines.o: $(SRCDIR)/conservationRoutines.h $(SRCDIR)/LP_ompi.h
$(OBJDIR)/EntropyCalculations.o: $(SRCDIR)/EntropyCalculations.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h
$(OBJDIR)/EquilibriumSolution.o: $(SRCDIR)/EquilibriumSolution.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/FieldCalculations.h #$(SRCDIR)/advection_1.h 
$(OBJDIR)/MarginalCreation.o: $(SRCDIR)/MarginalCreation.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h
$(OBJDIR)/MomentCalculations.o: $(SRCDIR)/MomentCalculations.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/FieldCalculations.h #$(SRCDIR)/advection_1.h
$(OBJDIR)/NegativityChecks.o: $(SRCDIR)/NegativityChecks.h  $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h
$(OBJDIR)/FieldCalculations.o: $(SRCDIR)/FieldCalculations.h $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h
$(OBJDIR)/SetInit_1.o: $(SRCDIR)/SetInit_1.h $(SRCDIR)/LP_ompi.h $(SRCDIR)/advection_1.h


wts: $(objects_wt)
	@echo "Building mpi weights"
	@$(MPICC) -O2 -openmp $(objects_wt) -o $(EXECDIR)weight.out  

clean:
	$(RM) $(OBJDIR)/*.o 
#	$(RM) $(EXECDIR)/*.out

# icpc -O2 -openmp LP_main.cpp -I$TACC_FFTW3_INC -I$TACC_MKL_INC -L$TACC_FFTW3_LIB -lfftw3_threads -lfftw3 -lpthread -lm -Wl,-rpath,$TACC_MKL_LIB -L$TACC_MKL_LIB -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread
