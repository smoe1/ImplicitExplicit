### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#
InitApp  =  InitApp
DFluxFunc = DFluxFunc

DogSolveUser  = $(DOGPACK)/lib/1d/DogSolveTDRK

# uncomment these lines to use the HLL(E) Riemann solver
#RiemannSolve     = $(LIB1D)/RiemannSolveHLL
#RiemannSolveLxW  = $(LIB1D)/LaxWendroff/RiemannSolveLxW_HLL

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
EulerParams.cpp \

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
