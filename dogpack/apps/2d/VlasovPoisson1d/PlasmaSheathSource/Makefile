
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

DogSolveUser     = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveSL

ComputeElecField = ComputeElecField
BeforeStep       = BeforeStep
ConSoln          = ConSoln
ApplyLimiter     = $(DOGPACK)/lib/2d/cart/ApplyLimiter_do_nothing

# Hard code zero boundary conditions
StepAdvec        = StepAdvec

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
