### definitions ###
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.defs

DogSolveUser           = $(DOGPACK)/apps/semi_lagrangian/lib/1d/DogSolveSL

# If we want periodic boundary conditions, why is there still a boundary
# conditions file sitting in this folder? (-DS)
SetBndValues           = SetBndValues_periodic

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES)

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
