### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#
InitApp  =  InitApp
DFluxFunc = DFluxFunc

DogSolveUser           = $(DOGPACK)/lib/1d/DogSolveTDRK

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
EulerParams.cpp \

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
