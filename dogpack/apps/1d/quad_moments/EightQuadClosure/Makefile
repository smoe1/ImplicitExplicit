### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#

InitApp      = InitApp
SetBndValues = SetBndValues_extrap


# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) \
QuadMomentParams.o \
GetMuts.o\
NewtonSys.o\
GetWeights.o\
LinTrans.o


ALL_SOURCES = $(COMMON_SOURCES) \
QuadMomentParams.cpp \
GetMuts.cpp\
NewtonSys.cpp\
GetWeights.cpp\
LinTrans.cpp


### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
