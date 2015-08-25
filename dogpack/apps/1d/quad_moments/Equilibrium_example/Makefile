### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#

include $(DOGPACK)/apps/1d/Poisson/lib/Makefile.defs

InitApp            = InitApp
BeforeQinit        = BeforeQinit
BeforeFullTimeStep = BeforeFullTimeStep
AfterFullTimeStep  = AfterFullTimeStep
SetBndValues       = SetBndValues_periodic
ApplyLimiter       = ApplyLimiter
ConSoln            = ConSoln
RiemannSolve       = RiemannSolveKinetic

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(POISSON1D_LIB_OBJECTS) \
QuadMomentParams.o \
MomentInversion.o \
GetAlpha.o \
ComputeElecField.o 

ALL_SOURCES = $(COMMON_SOURCES) $(POISSON1D_LIB_SOURCES) \
QuadMomentParams.cpp \
MomentInversion.cpp \
GetAlpha.cpp \
ComputeElecField.cpp

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
