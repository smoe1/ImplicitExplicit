### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# compiler options
#
#FFLAGS =  -c -O4
#LFLAGS =  
# -show: display g++ commands that would be executed
# -echo: show *all* commands and execute them
#LINK    = /usr/local/hdf5/bin/h5c++ -show # -E

include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

# place to add files
#
ALL_OBJECTS = $(SL_LIB_OBJECTS)
ALL_SOURCES = $(SL_LIB_SOURCES) $(LIBSOURCES_COMMON)

### targets ###

all: \
  $(ALL_OBJECTS)

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
