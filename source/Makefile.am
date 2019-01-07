bin_PROGRAMS  = solver
AM_CPPFLAGS   = $(FFTW_CFLAGS) 
AM_CPPFLAGS  += -I$(OPENBLAS_INC)
LIBS          = $(BLAS_LIBS) $(MKL_LIBS)

h_sources   = EntropyCalculations.h LP_ompi.h NegativityChecks.h collisionRoutines_1.h \
	      EquilibriumSolution.h MarginalCreation.h SetInit_1.h conservationRoutines.h \
	      FieldCalculations.h MomentCalculations.h advection_1.h InputParsing.h

cpp_sources = EntropyCalculations.cpp LP_ompi.cpp NegativityChecks.cpp  collisionRoutines_1.cpp \
              EquilibriumSolution.cpp MarginalCreation.cpp SetInit_1.cpp conservationRoutines.cpp \
	      FieldCalculations.cpp MomentCalculations.cpp advection_1.cpp InputParsing.cpp

solver_SOURCES = $(cpp_sources) $(h_sources)
