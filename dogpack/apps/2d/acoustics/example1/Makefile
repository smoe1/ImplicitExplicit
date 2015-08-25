
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# overridden library objects
#
# InitApp = InitApp

# place to add object modules
#
APP_OBJECTS = \
	$(DOGPACK)/apps/2d/acoustics/example1/AcousticParams.o \
	$(DOGPACK)/apps/2d/acoustics/example1/AppSolver.o \

APP_SOURCES = \
	$(DOGPACK)/apps/2d/acoustics/example1/AcousticParams.cpp \
	$(DOGPACK)/apps/2d/acoustics/example1/AppSolver.cpp \
	*.h \

# Two-derivative solver:
DogSolveUser        = $(DOGPACK)/lib/2d/cart/DogSolveTD

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##


# recombobulate stuff ...
ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(APP_SOURCES)

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
