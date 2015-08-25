### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# compiler options
#
include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.defs

# place to add files
#
ALL_OBJECTS = $(SL_LIB_OBJECTS)
ALL_SOURCES = $(SL_LIB_SOURCES) $(LIBSOURCES_COMMON)

### targets ###

lib: \
  $(ALL_OBJECTS)

#include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.targets
include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
