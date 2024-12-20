#===============================================================
# Project Settings
#===============================================================

PROJECT_NAME = "FEMaster"

#===============================================================
# Compiler Options and Configurations
#===============================================================

# Compiler settings
NVCC = nvcc
CXX  = g++


#===============================================================
# Default Feature Flags
#===============================================================

openmp           ?= 1   # Enable/disable OpenMP
mkl              ?= 0   # Enable/disable MKL (Math Kernel Library)
mkl_sequential   ?= 0   # Enable/disable Sequential MKL
cuda_dp          ?= 1   # Enable/disable CUDA Double Precision
debug            ?= 0   # Enable/disable Debug mode
double_precision ?= 1 # Enable/disable Double Precision
#===============================================================
# General Compiler Flags
#===============================================================

WARNFLAGS  =
CXXFLAGS   = -std=c++17 -O3 -I /usr/local/include $(WARNFLAGS)
NVCCFLAGS  = -std=c++17 -O3 -I /usr/local/include --expt-relaxed-constexpr $(WARNFLAGS)

#===============================================================
# Custom Feature Flags (Conditional Flags)
#===============================================================

# OpenMP support
CXXFLAGS  += $(if $(filter 1,$(openmp)),-fopenmp)
NVCCFLAGS += $(if $(filter 1,$(openmp)),-Xcompiler=-fopenmp)

# Ensure MKL is enabled when sequential MKL is requested
mkl := $(if $(filter 1,$(mkl_sequential)),1,$(mkl))

# MKL support
FEATURE_FLAGS += $(if $(filter 1,$(mkl)),-DMKL_LP64)
FEATURE_FLAGS += $(if $(filter 1,$(mkl)),-DUSE_MKL)

# Sequential or parallel MKL
FEATURE_FLAGS += $(if $(filter 1,$(mkl_sequential)),-DUSE_MKL_SEQUENTIAL)
LIBS          += $(if $(filter 1,$(mkl_sequential)),$(MKL_LIBS),$(if $(filter 1,$(mkl)),$(MKL_LIBS)))

# Debug mode
FEATURE_FLAGS += $(if $(filter 1,$(debug)),-DNDEBUG)
CXXFLAGS      += $(if $(filter 1,$(debug)),-g)
NVCCFLAGS     += $(if $(filter 1,$(debug)),-G -g)

# cuda double precision
FEATURE_FLAGS += $(if $(filter 1,$(cuda_dp)),-DCUDA_DOUBLE_PRECISION)

# double precision
FEATURE_FLAGS += $(if $(filter 1,$(double_precision)),-DDOUBLE_PRECISION)


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

# Change the executable definitions to remove .exe
EXE_CPU  := $(BINDIR)/$(PROJECT_NAME)
EXE_GPU  := $(BINDIR)/$(PROJECT_NAME)_gpu
EXE_TST  := $(BINDIR)/$(PROJECT_NAME)_test

# Add these new utility targets at the end of your makefile
info:
	@echo "Build Configuration:"
	@echo "  Project        : $(PROJECT_NAME)"
	@echo "  Platform       : $(UNAME) ($(shell uname -m))"
	@echo "  Compiler       : $(COMPILER_VERSION)"
	@echo "  C++ Flags      : $(CXXFLAGS)"
	@echo "  NVCC Flags     : $(NVCCFLAGS)"
	@echo "  MKL Enabled    : $(mkl)"
	@echo "  OpenMP Enabled : $(openmp)"
	@echo "  Debug Mode     : $(debug)"

clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(OBJDIR) $(BINDIR)

help:
	@echo "Available targets:"
	@echo "  all    : Build everything (default)"
	@echo "  cpu    : Build CPU-only version"
	@echo "  gpu    : Build GPU-enabled version"
	@echo "  tests  : Build test suite"
	@echo "  clean  : Remove build artifacts"
	@echo "  info   : Show build configuration"
	@echo "  help   : Show this help message"

#===============================================================
# Build Targets
#===============================================================

all: info cpu gpu tests

cpu: info $(EXE_CPU)

gpu: CXXFLAGS  += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: info $(EXE_GPU) clean-exp-lib

tests: info $(EXE_TST)

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
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) $^ $(TEST_LIBS) -o $@

# Object generation rules with timing
$(GPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s.%N); \
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s.%N); \
	DURATION=$$(echo "scale=1; $$END - $$START" | bc); \
	printf "%-50s : %6.1f seconds\n" "Compiling $<" $$DURATION

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s.%N); \
	$(NVCC) $(NVCCFLAGS) $(FEATURE_FLAGS) -x cu -c $< -o $@; \
	END=$$(date +%s.%N); \
	DURATION=$$(echo "scale=1; $$END - $$START" | bc); \
	printf "%-50s : %6.1f seconds\n" "Compiling $<" $$DURATION

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	@START=$$(date +%s.%N); \
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s.%N); \
	DURATION=$$(echo "scale=1; $$END - $$START" | bc); \
	printf "%-50s : %6.1f seconds\n" "Compiling $<" $$DURATION

$(GPP_OBJDIR)/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s.%N); \
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s.%N); \
	DURATION=$$(echo "scale=1; $$END - $$START" | bc); \
	printf "%-50s : %6.1f seconds\n" "Compiling $<" $$DURATION


# Include the generated dependency files
-include $(GPP_OBJS:.o=.d)
-include $(GPU_CPP_OBJS:.o=.d)
-include $(CU_OBJS:.o=.d)

#===============================================================
# Clean Targets
#===============================================================

clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib

