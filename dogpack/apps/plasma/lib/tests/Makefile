
include $(DOGPACK)/lib/Makefile.defs # inherit interface
include $(DOGPACK)/apps/plasma/lib/Makefile.defs
include $(DOGPACK)/apps/plasma/lib/sys/tests/Makefile.defs

Source2Fluid_OBJECTS = \
  $(POLYNOMIAL_OBJECTS) \
  $(PLASMA_LIB)/Source2Fluid.o \
  $(LIB_SYS)/assert.h \
  $(LIB_SYS)/assert.cpp \
  test_Source2Fluid.o \

Source2Fluid_SOURCES = \
  $(POLYNOMIAL_SOURCES) \
  $(PLASMA_LIB)/Source2Fluid.h \
  $(PLASMA_LIB)/Source2Fluid.cpp \
  test_Source2Fluid.cpp \

all: test_Source2Fluid

test_Source2Fluid: $(Source2Fluid_OBJECTS) $(Source2Fluid_SOURCES) 
	$(LINK) $(Source2Fluid_SOURCES) -o test_Source2Fluid

clean: clean_Source2Fluid

clean_Source2Fluid:
	rm -f $(Source2Fluid_OBJECTS) test_Source2Fluid

include $(DOGPACK)/lib/sys/Makefile.targets
