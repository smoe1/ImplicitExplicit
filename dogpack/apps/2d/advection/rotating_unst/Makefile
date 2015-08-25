
### definitions ###

# default/library definitions
COMPILE_WITH_UNST = 1
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
SetBndValues_Unst = SetBndValues_Unst
AfterQinit_Unst   = AfterQinit_Unst

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##


ALL_OBJECTS = $(COMMON_OBJECTS) NonConsFunc.o
ALL_SOURCES = $(COMMON_SOURCES) NonConsFunc.cpp

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
