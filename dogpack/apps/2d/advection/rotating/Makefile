
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

## ------------------------------------------------------------------------- ##
## Choose which DogSolveUser method to try:
##
##   Options include: semi-lagrangian scheme, hybrid semi-lagrangian, and
##   two-derivative method.
##
#DogSolveUser       = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveSL
#DogSolveUser       = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveHybrid
DogSolveUser        = $(DOGPACK)/lib/2d/cart/DogSolveTD
## ------------------------------------------------------------------------- ##

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
#D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##


ALL_OBJECTS = \
$(COMMON_OBJECTS) \
$(SL_LIB_OBJECTS) \
AppSolver.o \
NonConsFunc.o \

ALL_SOURCES = \
$(COMMON_SOURCES) \
$(SL_LIB_SOURCES) \
AppSolver.cpp \
NonConsFunc.cpp \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
