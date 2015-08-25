
include $(DOGPACK)/lib/sys/Makefile.defs # inherit interface
include $(DOGPACK)/apps/plasma/lib/sys/Makefile.defs
include $(DOGPACK)/apps/plasma/lib/sys/tests/Makefile.defs

ALL_SOURCES = \
  $(POLYNOMIAL_SOURCES) \

ALL_OBJECTS = \
  $(POLYNOMIAL_OBJECTS) \

CLEAN_TARGETS = clean_Polynomial

test_Polynomial:  $(POLYNOMIAL_OBJECTS) $(POLYNOMIAL_SOURCES)
	$(LINK) $(POLYNOMIAL_OBJECTS) -o test_Polynomial

all: test_Polynomial

clean_Polynomial:
	rm -f $(POLYNOMIAL_OBJECTS) test_Polynomial

include $(DOGPACK)/lib/sys/Makefile.targets

