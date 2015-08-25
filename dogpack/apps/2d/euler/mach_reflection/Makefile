### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
AfterQinit   = AfterQinit
RiemannSolve = $(DOGPACK)/lib/2d/RiemannSolveLLF

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:

DogSolveUser        = $(DOGPACK)/lib/2d/cart/DogSolveTD
DFluxFunc           = DFluxFunc

## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##


ALL_OBJECTS = $(COMMON_OBJECTS) \
AppSolver.o \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
AppSolver.cpp \
EulerParams.cpp \
*.h \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
