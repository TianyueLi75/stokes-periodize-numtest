# Example Makefile for building CSBQ projects

PVFMM_DIR ?= ./extern/pvfmm
-include $(PVFMM_DIR)/MakeVariables

# Directories for SCTL includes and quadrature tables
SCTL_INCLUDE_DIR ?= $(PVFMM_DIR)/SCTL/include 
CSBQ_INCLUDE_DIR ?= ./extern/CSBQ/include
SCTL_DATA_PATH ?= ./extern/CSBQ/data

# Compiler settings
CXX = $(CXX_PVFMM) 
CXXFLAGS = -std=c++17 -fopenmp # Need C++11 and OpenMP

# Define the path for quadrature tables and enable quadruple precision (for reading quadrature tables)
CXXFLAGS += -DSCTL_DATA_PATH=$(SCTL_DATA_PATH)
CXXFLAGS += -DSCTL_QUAD_T=__float128


############################################
#             OPTIONAL FLAGS               #
############################################

# Optional flags for debugging or release
DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CXXFLAGS += -O2 -g -fno-omit-frame-pointer -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract,float-divide-by-zero,float-cast-overflow \
                -fno-sanitize-recover=all -fstack-protector # Debug build
    CXXFLAGS += -DSCTL_MEMDEBUG # Enable SCTL memory checks
else
    CXXFLAGS += -O3 -march=native -DNDEBUG
endif

# Enable warnings
CXXFLAGS += -Wall -Wfloat-conversion

# Enable profiling
CXXFLAGS += -DSCTL_PROFILE=100 -DSCTL_VERBOSE
# CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE

# Enable MPI (CXX must be set to mpicxx)
CXXFLAGS += -DSCTL_HAVE_MPI

# Enable BLAS and LAPACK
# CXXFLAGS += -lblas -DSCTL_HAVE_BLAS                        # Use BLAS
# CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK                    # Use LAPACK
# CXXFLAGS += -lopenblas -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # Use OpenBLAS
# CXXFLAGS += -mkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK       # Use MKL BLAS and LAPACK (Intel compiler)
# CXXFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # Use MKL BLAS and LAPACK (non-Intel compiler)
# CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK
CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK

# Enable FFTW
# CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
# CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
# CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL

# Enable SVML (Intel compiler)
# CXXFLAGS += -DSCTL_HAVE_SVML

# Enable PVFMM
PVFMM_INC_DIR = $(PVFMM_DIR)/include
PVFMM_LIB_DIR = $(PVFMM_DIR)/lib/.libs
CXXFLAGS += -DSCTL_HAVE_PVFMM -I$(PVFMM_INC_DIR)
LDLIBS += $(PVFMM_LIB_DIR)/libpvfmm.a

# Settings for stack trace
ifeq ($(shell uname -s),Darwin)
    CXXFLAGS += -g -rdynamic -Wl,-no_pie # For stack trace on Mac
else
    CXXFLAGS += -gdwarf-4 -g -rdynamic # For stack trace on other systems
endif

# Utility commands
RM = rm -f
MKDIRS = mkdir -p

# Directories for binaries, objects, includes, tests
BINDIR = ./bin
OBJDIR = ./obj
INCDIR = ./include
TESTDIR = ./test

TEST_BIN = \
    $(BINDIR)/examples \
    $(BINDIR)/periodization_time \
	$(BINDIR)/precompute_time \
	$(BINDIR)/test_convergence \
	$(BINDIR)/test_periodicity \
	$(BINDIR)/test2_ptcl_conv \
	$(BINDIR)/test2 \

# Test target: build all test binaries
test: $(TEST_BIN)

# Example script with visualizations in singly- and doubly-periodicity
examples: $(BINDIR)/examples

# Timing call for scaling; also includes visualization.
timing: ${BINDIR}/test2

# Timing call for precomputing periodizing operators
timing_precomp: ${BINDIR}/precompute_time

# Timing call for checking periodization overhead
timing_periodization: ${BINDIR}/periodization_time

# Test accuracy of singly-periodic code through manufactured solutions on singly-periodic spherical suspensions
test_manufactured_soln: $(BINDIR)/test2_ptcl_conv
# Test self convergence of solver in singly, doubly, or triply periodic geometries
test_selfconv: $(BINDIR)/test_convergence

# Verification: check periodicity at two sides of periodic box.
test_peri: $(BINDIR)/test_periodicity



# Rules for building binaries
$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $^ $(LDLIBS) -o $@ $(CXXFLAGS)
ifeq "$(OS)" "Darwin"
	/usr/bin/dsymutil $@ -o $@.dSYM
endif

# Rules for compiling test source files
$(OBJDIR)/%.o: $(TESTDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(SCTL_INCLUDE_DIR) -I$(CSBQ_INCLUDE_DIR) -c $^ -o $@

# Clean target: remove all binaries and object files
clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*

