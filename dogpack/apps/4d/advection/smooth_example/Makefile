### definitions ###

# default/library definitions
include $(DOGPACK)/lib/4d/Makefile.defs

# compiler options -- see file $(DOGPACK)/config/Makefile.config

# overridden library objects
SetBndValues  = SetBndValues
#ConstructL = ConstructL

ALL_OBJECTS = \
$(COMMON_OBJECTS) \
AppSolver.o \
InitialParams.o \

ALL_SOURCES = \
$(COMMON_SOURCES) \
AppSolver.cpp \
InitialParams.cpp \

### targets ###

include $(DOGPACK)/lib/4d/Makefile.targets

### DO NOT remove this line - make depends on it ###
