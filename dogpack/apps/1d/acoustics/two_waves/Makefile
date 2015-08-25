### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# Semi-lagrangian library, needed to run the semi-Lagrangian option in the
# code.
include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.defs

# overridden library objects
#
#DogSolveUser           = $(DOGPACK)/lib/1d/DogSolveTDRK
DogSolveUser           = DogSolveSL

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES)

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
