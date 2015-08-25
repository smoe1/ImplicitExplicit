### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#
InitApp     =  InitApp
AfterQinit  =  AfterQinit

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) \
MHDParams.o \
InitialParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
MHDParams.cpp \
InitialParams.cpp \

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
