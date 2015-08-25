
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# compiler options -- see file $(DOGPACK)/config/Makefile.config

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

### Choose Which Solving method to use if desired ###
#DogSolveUser       = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveSL
DogSolveUser       = $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/DogSolveHybrid

HybridSourceTermFunc1D = HybridSourceTermFunc1D

#IntegratePhi       = IntegratePhi
#InitSLState        = InitSLState

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
