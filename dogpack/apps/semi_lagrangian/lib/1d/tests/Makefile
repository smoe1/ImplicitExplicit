include $(DOGPACK)/lib/Makefile.defs
include $(DOGPACK)/apps/semi_lagrangian/lib/1d/Makefile.defs

DOGDEFS_SOURCES = \
 $(TOPLIB)/constants.h \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(TOPLIB)/tensors.h \
 $(TOPLIB)/tensors.cpp \
 $(DOGPACK)/apps/semi_lagrangian/lib/1d/dg_cell.h \
 $(DOGPACK)/apps/semi_lagrangian/lib/1d/dg_cell.cpp \

DOGDEFS_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/debug.o \
 $(TOPLIB)/tensors.o \
 $(DOGPACK)/apps/semi_lagrangian/lib/1d/dg_cell.o \

DOGDEFS_DEPENDENCIES = \
 $(DOGDEFS_SOURCES) \
 $(DOGDEFS_OBJECTS) \

TEST_DG_CELL_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_dg_cell.o

%.o : %.cpp
	$(LINK) $(FFLAGS) $< -o $@

EXECUTABLES = \
  test_dg_cell

all: $(EXECUTABLES)

test_dg_cell:  $(DOGDEFS_DEPENDENCIES) test_dg_cell.cpp test_dg_cell.o
	$(LINK) $(LFLAGS) $(TEST_DG_CELL_OBJECTS) -o test_dg_cell

CLEAN_FILES = *.o $(EXECUTABLES)

include $(DOGPACK)/lib/Makefile.targets
