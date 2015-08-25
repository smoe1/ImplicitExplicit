
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

DogSolveUser       = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveSL
SetAdvecSpeed      = SetAdvecSpeed
SetAdvecSpeedExtra = SetAdvecSpeedExtra
InitSLState        = InitSLState

ALL_OBJECTS = \
$(COMMON_OBJECTS) \
$(SL_LIB_OBJECTS) \
AppSolver.o \

ALL_SOURCES = \
$(COMMON_SOURCES) \
$(SL_LIB_SOURCES) \
AppSolver.cpp \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
