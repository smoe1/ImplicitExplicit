
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# place to add object modules
#
InitApp         = InitApp
AuxFunc         = $(DOGPACK)/lib/1d/blanks/AuxFunc
FluxFunc        = $(DOGPACK)/lib/1d/blanks/FluxFunc
QinitFunc       = $(DOGPACK)/lib/1d/blanks/QinitFunc
ProjectLeftEig  = $(DOGPACK)/lib/1d/blanks/ProjectLeftEig
ProjectRightEig = $(DOGPACK)/lib/1d/blanks/ProjectRightEig
SetBndValues    = $(DOGPACK)/lib/1d/blanks/SetBndValues
SetWaveSpd      = $(DOGPACK)/lib/1d/blanks/SetWaveSpd

include $(DOGPACK)/apps/1d/Poisson/lib/Makefile.defs

ALL_OBJECTS = $(COMMON_OBJECTS) $(POISSON1D_LIB_OBJECTS)

ALL_SOURCES = $(COMMON_SOURCES) $(POISSON1D_LIB_SOURCES)

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
