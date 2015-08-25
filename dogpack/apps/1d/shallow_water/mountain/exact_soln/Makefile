
### definitions ###

# default/library definitions
include $(DOGPACK)/lib/1d/Makefile.defs

### Substitutions ###
AuxFunc         = $(DOGPACK)/lib/1d/blanks/AuxFunc
FluxFunc        = $(DOGPACK)/lib/1d/blanks/FluxFunc
QinitFunc       = $(DOGPACK)/lib/1d/blanks/QinitFunc
ProjectLeftEig  = $(DOGPACK)/lib/1d/blanks/ProjectLeftEig
ProjectRightEig = $(DOGPACK)/lib/1d/blanks/ProjectRightEig
SetBndValues    = $(DOGPACK)/lib/1d/blanks/SetBndValues
SetWaveSpd      = $(DOGPACK)/lib/1d/blanks/SetWaveSpd
SourceTermFunc  = $(DOGPACK)/lib/1d/blanks/SourceTermFunc

### targets ###

shllw.exe: $(ALL_OBJECTS)
	$(LINK) $(ALL_OBJECTS) -o shllw.exe

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
