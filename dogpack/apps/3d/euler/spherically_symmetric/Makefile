
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/3d/Makefile.defs

# compiler options -- see file $(DOGPACK)/config/Makefile.config

# overridden library objects
SetBndValues         = SetBndValues
AfterQinit           = AfterQinit
L2ProjectInitialCond = L2ProjectInitialCond

ALL_OBJECTS = \
$(COMMON_OBJECTS) \
AppSolver.o \
EulerParams.o \
InitialParams.o \

ALL_SOURCES = \
$(COMMON_SOURCES) \
AppSolver.cpp \
EulerParams.cpp \
InitialParams.cpp \

### targets ###

include $(DOGPACK)/lib/3d/Makefile.targets

### DO NOT remove this line - make depends on it ###
