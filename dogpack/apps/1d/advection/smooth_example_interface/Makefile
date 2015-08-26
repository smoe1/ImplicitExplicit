### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs


# compiler options
#
# overridden library objects
#

# Semi-lagrangian library, needed to run the semi-Lagrangian option in the
# code.
include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.defs

## ------------------------------------------------------------------------- ##
# List of a few different options explored for this problem:
#
#   DogSolveRK4  : classical, fourth-order time integration.
#   DogSolveSL   : Semi-Lagrangian solver
#   DogSolveTDRK : Fourth-order, two-derivative solver
#
# We currently have no mechanism in place to select out of these various
# solvers, so this Makefile has to link to the appropriate solver.
#

AuxFunc     = AuxFunc
BeforeFullTimeStep = BeforeFullTimeStep

#DogSolveUser           = $(DOGPACK)/apps/semi_lagrangian/lib/1d/DogSolveSL
#DogSolveUser           = $(DOGPACK)/lib/1d/DogSolveRK4
#AfterQinit             = AfterQinit
ApplyLimiter           =ApplyLimiter
#RiemannSolve          = $(DOGPACK)/lib/1d/RiemannSolveHLL

## ------------------------------------------------------------------------- ##

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS) \
GetCFL_Interface.o \

ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES) \
GetCFL_Interface.cpp \

### targets ###

include $(DOGPACK)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
