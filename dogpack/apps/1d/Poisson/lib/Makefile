### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# compiler options
#
#FFLAGS =  -c -O4
#LFLAGS =  
# -show: display g++ commands that would be executed
# -echo: show *all* commands and execute them
#LINK    = /usr/local/hdf5/bin/h5c++ -show # -E

include $(DOGPACK)/apps/1d/Poisson/lib/Makefile.defs

# place to add files
#
ALL_OBJECTS = $(POISSON1D_LIB_OBJECTS)
ALL_SOURCES = $(POISSON1D_LIB_SOURCES) $(LIBSOURCES_COMMON)

### targets ###

poisson1dlib: \
  $(POISSON1D_LIB_OBJECTS)

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
