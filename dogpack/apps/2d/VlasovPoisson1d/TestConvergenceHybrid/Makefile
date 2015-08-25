
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# compiler options (See $(DOGPACK)/config/Makefile.config)

VLASOV_LIB  = -I$(DOGPACK)/apps/2d/VlassovPoisson1d/lib
POISSON_LIB = -I$(DOGPACK)/apps/1d/Poisson/

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs
include $(DOGPACK)/apps/2d/VlasovPoisson1d/lib/Makefile.defs

### Do not change this to DogSolveSL.  See TestConvergenceHybrid for that problem
DogSolveUser     = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveHybrid

BeforeStep       = BeforeStep
AfterStep        = AfterStep

### These routines are particular to the TestConvergence Problem: ###
SetAdvecSpeedExtra     = SetAdvecSpeedExtra
InitHybridState        = InitHybridState
IntegratePsi           = IntegratePsi
HybridSourceTermFunc1D = HybridSourceTermFunc1D

ALL_OBJECTS = $(COMMON_OBJECTS) \
$(SL_LIB_OBJECTS) \
$(VLASOV_LIB_OBJECTS) \
AppSolver.o \
#
ALL_SOURCES = $(COMMON_SOURCES) \
$(SL_LIB_SOURCES) \
$(VLASOV_LIB_SOURCES) \
AppSolver.cpp \
#
### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
