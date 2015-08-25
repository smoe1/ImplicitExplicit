
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
AfterQinit  =  AfterQinit

# Lax-Wendroff methods
DFluxFunc       = DFluxFunc
D2FluxFunc      = D2FluxFunc

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
