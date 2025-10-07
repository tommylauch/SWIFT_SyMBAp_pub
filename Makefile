# Makefile for SWIFT N-body integration package

# Configuration
FC = gfortran
FFLAGS = -O2 -ftree-vectorize -march=native
CPP = cpp
CPPFLAGS = -U_FXDR_AVAIL

# change these if needed
USE_OPENMP = 1

# Directories and files
LIB = libswift.a
STANDARD_DIRS = anal bs coord discard io lyap lyap2 orbel rmvs rmvs2 rmvs3 rmvs4 tu4 obl util helio symba5 skeel \
                mvs/drift mvs/getacch mvs/kickvh mvs/step mvs_array/drift mvs_array/getacch mvs_array/kickvh

# Programs (with full paths)
ifeq ($(USE_OPENMP),1)
	OMP_DIRS = symba5p
	MAIN = $(patsubst %.f,%,$(wildcard main/swift_*.f))
else
	OMP_DIRS = 
	MAIN = $(filter-out main/swift_symba5p,$(patsubst %.f,%,$(wildcard main/swift_*.f)))
endif
TOOL = $(patsubst %.f,%,$(wildcard tools/*.f))

# Source files and objects
SRC = $(foreach dir,$(STANDARD_DIRS) $(OMP_DIRS),$(wildcard $(dir)/*.f) $(wildcard $(dir)/*.F))
OBJS = $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(SRC)))

# Include files
INC_FILES = swift.inc symba5/symba5.inc symba5p/symba5p.inc

# Build rules
.PHONY: main tool library clean help
.SUFFIXES: .f .F .o

main: $(SRC) $(INC_FILES)
	@$(MAKE) --no-print-directory $(MAIN)
	@echo "make done"

tool: $(SRC) $(INC_FILES)
	@$(MAKE) --no-print-directory $(TOOL)
	@echo "make done"

library: $(SRC) $(INC_FILES)
	@$(MAKE) --no-print-directory $(LIB)
	@echo "make done"

$(LIB): $(SRC) $(INC_FILES)
	@rm -f $@
	@$(MAKE) -s --no-print-directory $(OBJS)
	@ar vr $@ $(OBJS) > /dev/null 2>&1
	@ranlib $@
	@rm -f $(OBJS)
	@find . -name "*_CPP.f" -delete
	@echo "$(LIB) built successfully"

# Library compilation
%.o: %.f $(INC_FILES)
	@$(FC) $(FFLAGS) -I./util -c $< -o $@
# 	Preprocessing compilation
%.o: %.F $(INC_FILES)
	@$(CPP) $(CPPFLAGS) $< $*_CPP.f
	@$(FC) $(FFLAGS) -I./util -c $*_CPP.f -o $@
	@rm -f $*_CPP.f
# 	Special rule for util_signal.o to keep .mod file in util/ directory
util/util_signal.o: util/util_signal.f $(INC_FILES)
	@cd util && $(FC) $(FFLAGS) -I. -c util_signal.f -o util_signal.o
# 	Special rule for OpenMP files
symba5p/%.o: symba5p/%.f $(INC_FILES)
	@$(FC) $(FFLAGS) -I./util -fopenmp -c $< -o $@

# main program compilation
$(filter-out main/swift_symba5p,$(MAIN)): %: %.f $(LIB) $(INC_FILES)
	@echo "--- MAIN: Compiling $(notdir $*) ---"
	@$(FC) $(FFLAGS) -I./util -o $@ $< -L. -lswift
# 	Special rule for OpenMP files
main/swift_symba5p: main/swift_symba5p.f $(LIB) $(INC_FILES)
	@echo "--- MAIN: Compiling swift_symba5p ---"
	@$(FC) $(FFLAGS) -I./util -fopenmp -o $@ $< -L. -lswift

# tools program compilation
$(TOOL): %: %.f $(LIB) $(INC_FILES)
	@echo "--- TOOL: Compiling $(notdir $*) ---"
	@$(FC) $(FFLAGS) -I./util -o $@ $< -L. -lswift

# Clean
clean:
	@rm -f $(OBJS) $(LIB) ./util/util_signal.mod
	@find . -name "*_CPP.f" -delete
	@for prog in $(MAIN); do rm -f $$prog; done
	@for prog in $(TOOL); do rm -f $$prog; done

# Help
help:
	@echo "SWIFT N-body Integration Package - Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  main       - Build libswift.a and main (default)"
	@echo "  tools      - Build libswift.a and standalone tools in tools/"
	@echo "  library    - Build libswift.a only"
	@echo "  clean      - Remove all generated files"
	@echo ""
	@echo "Configuration:"
	@echo "  FC = $(FC)"
	@echo "  FFLAGS = $(FFLAGS)"
	@echo "  USE_OPENMP = $(USE_OPENMP) (1=enabled, 0=disabled)"
	@echo ""
	@echo "Build options:"
	@echo "  help       - Show this help message"
	@echo "  make -j    - Parallel build"
	@echo "  make USE_OPENMP=0 - Build without OpenMP support"
	@echo ""
	@echo "Main programs that will be built:"
	@for prog in $(MAIN); do echo "  $$(basename $$prog)"; done
	@echo ""
	@echo "Tools that will be built:"
	@for prog in $(TOOL); do echo "  $$(basename $$prog)"; done
