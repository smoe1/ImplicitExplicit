
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs

# compiler options -- see file $(DOGPACK)/config/Makefile.config

# overridden library objects
#
include $(DOGPACK)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

AfterQinit          = AfterQinit

### Choose Which Solving method to use if desired ###
DogSolveUser        = DogSolveUser

ALL_OBJECTS = \
$(COMMON_OBJECTS) \
AppSolver.o \
StepAdvecCC.o \

ALL_SOURCES = \
$(COMMON_SOURCES) \
AppSolver.cpp \
StepAdvecCC.cpp \

### targets ###

include $(DOGPACK)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
