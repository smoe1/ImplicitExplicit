### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

# overridden library objects
#
DogSolveUser  = $(DOGPACK)/lib/1d/DogSolveTDRK
#DogSolveUser  = $(DOGPACK)/lib/1d/DogSolveRK4
DFluxFunc     =  DFluxFunc
D2FluxFunc    =  D2FluxFunc

RiemannSolve = $(DOGPACK)/lib/1d/RiemannSolveLLF

ConSoln       = ConSoln

# place to add object modules
#
#ALL_OBJECTS = $(COMMON_OBJECTS) ...
#ALL_SOURCES = $(COMMON_SOURCES) ...

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
