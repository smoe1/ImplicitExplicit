
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/2d/Makefile.defs
include $(DOGPACK)/apps/plasma/lib/2d/Makefile.defs

# override library files with plasma library files
RiemannSolve       = $(PLASMA_LIB2D)/RiemannSolve
SetWaveSpd         = $(PLASMA_LIB2D)/SetWaveSpd

#LstarExtra = LstarExtra
#L2Project          = $(LIB2D_CART)/L2Project_divfree
#L2ProjectGrad      = $(LIB2D_CART)/L2ProjectGrad_divfree
ApplyLimiter       = $(PLASMA_LIB2D)/ApplyLimiter_divfree
#AfterFullTimeStep = AfterFullTimeStep
ConSoln = ConSoln
InitApp = InitApp

APP_OBJECTS = \
  App.o \

APP_SOURCES = \
  App.cpp \
  *.h \

ALL_OBJECTS = $(COMMON_OBJECTS) $(PLASMA_LIB2D_OBJECTS) $(APP_OBJECTS)
#
ALL_SOURCES = $(COMMON_SOURCES) $(PLASMA_LIB2D_SOURCES) $(APP_SOURCES)

### targets ###

include $(DOGPACK)/apps/plasma/lib/2d/Makefile.targets

problem: problem.tex
	pdflatex problem.tex

### DO NOT remove this line - make depends on it ###
