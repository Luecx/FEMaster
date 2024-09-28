#===============================================================
# Compiler Options and Configurations
#===============================================================

# Compiler settings
NVCC = nvcc
CXX  = g++

#===============================================================
# Default Feature Flags
#===============================================================

openmp         ?= 1   # Enable/disable OpenMP
mkl            ?= 0   # Enable/disable MKL (Math Kernel Library)
mkl_sequential ?= 0   # Enable/disable Sequential MKL
cuda_dp        ?= 1   # Enable/disable CUDA Double Precision
ar_pcs         ?= 0   # Enable/disable Show Array Processes
debug          ?= 0   # Enable/disable Debug mode

#===============================================================
# General Compiler Flags
#===============================================================

CXXFLAGS   = -std=c++17 -O3 -I /usr/local/include
NVCCFLAGS  = -std=c++17 -O3 -I /usr/local/include --expt-relaxed-constexpr
WARNFLAGS  = -Wall -Wno-maybe-uninitialized

#===============================================================
# Custom Feature Flags (Conditional Flags)
#===============================================================

# OpenMP support
ifeq ($(openmp), 1)
    CXXFLAGS  += -fopenmp
    NVCCFLAGS += -Xcompiler=-fopenmp
endif

# Ensure MKL is enabled when sequential MKL is requested
ifeq ($(mkl_sequential), 1)
    mkl = 1
endif

# MKL support
ifeq ($(mkl), 1)
    FEATURE_FLAGS += -DMKL_LP64
    FEATURE_FLAGS += -DUSE_MKL
    ifeq ($(mkl_sequential), 1)
        # Use sequential MKL
        FEATURE_FLAGS += -DUSE_MKL_SEQUENTIAL
        LIBS += $(MKL_LIBS)
    else
        # Use parallel MKL
        FEATURE_FLAGS += -DUSE_MKL_PARALLEL
        LIBS += $(MKL_LIBS)
    endif
endif

# CUDA Double Precision support
ifeq ($(cuda_dp), 1)
    FEATURE_FLAGS += -DCUDA_DOUBLE_PRECISION
endif

# Show Array Processes for debugging
ifeq ($(ar_pcs), 1)
    FEATURE_FLAGS += -DSHOW_ARRAY_PROCESSES
endif

# Debug mode
ifeq ($(debug), 1)
    FEATURE_FLAGS += -DNDEBUG
    CXXFLAGS += -g
    NVCCFLAGS += -G -g
endif

#===============================================================
# System Information (Optional Flags Based on System)
#===============================================================

UNAME := $(shell uname -s)

#===============================================================
# Test Libraries and Paths Based on System
#===============================================================

# Determine if the system is running on Apple Silicon
ifeq ($(UNAME),Darwin)
    ARCH := $(shell uname -m)
    ifeq ($(ARCH),arm64)
        TEST_LIBS = -L/opt/homebrew/lib -lgtest -lgtest_main -pthread
    else
        TEST_LIBS = -L/usr/local/lib -lgtest -lgtest_main -pthread
    endif
else
    TEST_LIBS = -L/usr/local/lib -lgtest -lgtest_main -pthread
endif

#===============================================================
# Compiler
#===============================================================

COMPILER_VERSION = $(shell $(CXX) --version | head -n 1)
# check if "clang" is inside the compiler version string
ifeq (,$(findstring clang,$(COMPILER_VERSION)))
	COMPILER = gcc
else
	COMPILER = clang
endif

#===============================================================
# MKL Library Paths and Libraries
#===============================================================

# Define MKL paths for different OS
MKL_PATH := $(MKLROOT)/lib/intel64
ifeq ($(OS),Darwin)
    MKL_PATH := $(MKLROOT)/lib
endif

# Define MKL library file (sequential or parallel)
MKL_LIB_FILE := libmkl_$(if $(filter 1,$(mkl_sequential)),sequential,intel_thread).a

# Compiler-specific MKL linking flags
MKL_COMMON_LIBS := $(MKL_PATH)/libmkl_intel_lp64.a $(MKL_PATH)/libmkl_core.a -lpthread -lm -ldl
ifeq ($(COMPILER), clang)
    MKL_LIBS := -Wl,-force_load,$(MKL_COMMON_LIBS) -Wl,-force_load,$(MKL_PATH)/$(MKL_LIB_FILE) $(if $(filter 0,$(mkl_sequential)), -L$(MKL_PATH) -liomp5)
else
    MKL_LIBS := -Wl,--start-group $(MKL_COMMON_LIBS) $(MKL_PATH)/$(MKL_LIB_FILE) -Wl,--end-group $(if $(filter 0,$(mkl_sequential)), -L$(MKL_PATH) -liomp5)
endif

#===============================================================
# CUDA Libraries
#===============================================================

NVCCLIBS := \
    -lcusolver \
    -lcublas \
    -lcusparse

#===============================================================
# Directories and File Paths
#===============================================================

SRCDIR      = src
TESTDIR     = tests
OBJDIR      = obj/
GPP_OBJDIR  = $(OBJDIR)gpp/
NVCC_OBJDIR = $(OBJDIR)nvcc/
CPP_OBJDIR  = $(NVCC_OBJDIR)cpp/
CU_OBJDIR   = $(NVCC_OBJDIR)cuda/
BINDIR      = bin

#===============================================================
# Source and Object Files
#===============================================================

# Define main source file and exclude it for test builds
MAIN_SRC    := $(SRCDIR)/main.cpp
CPP_SRCS    := $(filter-out $(MAIN_SRC), $(sort $(shell find $(SRCDIR) -name '*.cpp')))
CU_SRCS     := $(sort $(shell find $(SRCDIR) -name '*.cu'))
TST_SRCS    := $(sort $(shell find $(TESTDIR) -name '*.cpp'))

# Object files
GPP_OBJS     := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
GPU_CPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
CU_OBJS      := $(CU_SRCS:$(SRCDIR)/%.cu=$(CU_OBJDIR)/%.o)
TST_OBJS     := $(TST_SRCS:$(TESTDIR)/%.cpp=$(GPP_OBJDIR)/%.o) $(GPP_OBJS)

#===============================================================
# Executable Files
#===============================================================

EXE_CPU  := $(BINDIR)/FEMaster.exe
EXE_GPU  := $(BINDIR)/FEMaster_gpu.exe
EXE_TST  := $(BINDIR)/test_suite.exe

#===============================================================
# Display Compilation Flags
#===============================================================

show_flags:
	@echo "Compilation Flags:"
	@echo "-------------------"
	@echo "C++ Compiler : $(COMPILER_VERSION)"
	@echo "NVCC Compiler: $(NVCC)"
	@echo "CXXFLAGS     : $(CXXFLAGS)"
	@echo "NVCCFLAGS    : $(NVCCFLAGS)"
	@echo "LIBS         : $(LIBS)"
	@echo "Feature Flags: $(FEATURE_FLAGS)"
	@echo ""

#===============================================================
# Build Targets
#===============================================================

all: show_flags cpu gpu tests

cpu: show_flags $(EXE_CPU)

gpu: CXXFLAGS  += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: show_flags $(EXE_GPU) clean-exp-lib

tests: $(EXE_TST)

#===============================================================
# Build Rules
#===============================================================

CXXFLAGS   += -MMD -MP
NVCCFLAGS  += -MMD -MP

$(EXE_CPU): $(GPP_OBJS) $(MAIN_SRC:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) $^ $(LIBS) -o $@

$(EXE_GPU): $(GPU_CPP_OBJS) $(CU_OBJS) $(MAIN_SRC:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $^ $(LIBS) -o $@

$(EXE_TST): $(TST_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ $(TEST_LIBS) -o $@

# Object generation rules
$(GPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -c $< -o $@

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(FEATURE_FLAGS) -x cu -c $< -o $@

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $(FEATURE_FLAGS) -c $< -o $@

$(GPP_OBJDIR)/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Include the generated dependency files
-include $(GPP_OBJS:.o=.d)
-include $(GPU_CPP_OBJS:.o=.d)
-include $(CU_OBJS:.o=.d)

#===============================================================
# Clean Targets
#===============================================================

clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib

clean:
	rm -rf $(GPP_OBJDIR) $(NVCC_OBJDIR) $(BINDIR) $(OBJDIR)
