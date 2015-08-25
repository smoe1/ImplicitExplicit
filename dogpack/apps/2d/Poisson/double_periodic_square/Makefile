
### definitions ###

# default/library definitions
COMPILE_WITH_UNST = 1
include $(DOGPACK)/lib/2d/Makefile.defs
include $(DOGPACK)/apps/2d/Poisson/lib/unst/Makefile.defs

# locally redefined functions:
RhsFunc           = RhsFunc

ALL_OBJECTS = $(COMMON_OBJECTS) $(POISSON_2D_LIB_UNST_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(POISSON_2D_LIB_UNST_SOURCES)

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
