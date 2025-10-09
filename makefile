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

openmp             ?= 1   # Enable/disable OpenMP
mkl                ?= 0   # Enable/disable MKL (Math Kernel Library)
mkl_sequential     ?= 0   # Enable/disable Sequential MKL
cuda_dp            ?= 1   # Enable/disable CUDA Double Precision
debug              ?= 0   # Enable/disable Debug mode
double_precision   ?= 1   # Enable/disable Double Precision
eigen_fast_compile ?= 1   # Enable/disable Eigen fast compile mode
time_report        ?= 0   # Enable/disable time report for compilation

#===============================================================
# General Compiler Flags
#===============================================================

WARNFLAGS  =
CXXFLAGS   = -std=c++17 -O3 -I /usr/local/include $(WARNFLAGS)
NVCCFLAGS  = -std=c++17 -O3 -I /usr/local/include --expt-relaxed-constexpr $(WARNFLAGS)

ifneq (,$(shell which ccache 2>/dev/null))
  CXX  := ccache $(CXX)
  NVCC := ccache $(NVCC)
endif

#===============================================================
# Custom Feature Flags (Conditional Flags)
#===============================================================

# OpenMP support
CXXFLAGS  += $(if $(filter 1,$(openmp)),-fopenmp)
NVCCFLAGS += $(if $(filter 1,$(openmp)),-Xcompiler=-fopenmp) -diag-suppress 20014ws

# Ensure MKL is enabled when sequential MKL is requested
mkl := $(if $(filter 1,$(mkl_sequential)),1,$(mkl))

# MKL support
FEATURE_FLAGS += $(if $(filter 1,$(mkl)),-DMKL_LP64)
FEATURE_FLAGS += $(if $(filter 1,$(mkl)),-DUSE_MKL)

# Sequential or parallel MKL
FEATURE_FLAGS += $(if $(filter 1,$(mkl_sequential)),-DUSE_MKL_SEQUENTIAL)
LIBS          += $(if $(filter 1,$(mkl_sequential)),$(MKL_LIBS),$(if $(filter 1,$(mkl)),$(MKL_LIBS)))

# Debug mode
FEATURE_FLAGS += $(if $(filter 0,$(debug)),-DNDEBUG -DEIGEN_NO_DEBUG)
CXXFLAGS      += $(if $(filter 1,$(debug)),-g3,-g0)
NVCCFLAGS     += $(if $(filter 1,$(debug)),-G -g)

CXXFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(CXXFLAGS))
CXXFLAGS += $(if $(filter 1,$(debug)),-O2,-O3)
NVCCFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(NVCCFLAGS))
NVCCFLAGS += $(if $(filter 1,$(debug)),-O2,-O3)

# cuda double precision
FEATURE_FLAGS += $(if $(filter 1,$(cuda_dp)),-DCUDA_DOUBLE_PRECISION)

# double precision
FEATURE_FLAGS += $(if $(filter 1,$(double_precision)),-DDOUBLE_PRECISION)

ifeq ($(eigen_fast_compile),1)
  FEATURE_FLAGS += -DEIGEN_DONT_PARALLELIZE    # Threading via OpenMP/MKL statt Eigen
  FEATURE_FLAGS += -DEIGEN_UNROLLING_LIMIT=50  # weniger Aggro-Unrolling
  FEATURE_FLAGS += -DEIGEN_DONT_INLINE         # weniger Inlining -> schnelleres Kompilieren
endif

CXXFLAGS += $(if $(filter 1,$(time_report)),-ftime-report)

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
	@START=$$(date +%s%N); \
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s%N); \
	DURATION=$$(echo "scale=3; ($$END - $$START) / 1000000000" | bc); \
	printf "%-50s : %6s seconds\n" "Compiling $<" $$(echo $$DURATION | awk '{printf "%.1f", $$1}')

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(NVCC) $(NVCCFLAGS) $(FEATURE_FLAGS) -x cu -c $< -o $@; \
	END=$$(date +%s%N); \
	DURATION=$$(echo "scale=3; ($$END - $$START) / 1000000000" | bc); \
	printf "%-50s : %6s seconds\n" "Compiling $<" $$(echo $$DURATION | awk '{printf "%.1f", $$1}')

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s%N); \
	DURATION=$$(echo "scale=3; ($$END - $$START) / 1000000000" | bc); \
	printf "%-50s : %6s seconds\n" "Compiling $<" $$(echo $$DURATION | awk '{printf "%.1f", $$1}')

$(GPP_OBJDIR)/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -c $< -o $@; \
	END=$$(date +%s%N); \
	DURATION=$$(echo "scale=3; ($$END - $$START) / 1000000000" | bc); \
	printf "%-50s : %6s seconds\n" "Compiling $<" $$(echo $$DURATION | awk '{printf "%.1f", $$1}')


# Include the generated dependency files
-include $(GPP_OBJS:.o=.d)
-include $(GPU_CPP_OBJS:.o=.d)
-include $(CU_OBJS:.o=.d)


#===============================================================
# IWYU (Include-What-You-Use) Integration  [CPU-only version]
#===============================================================


# Detect IWYU tools
IWYU           ?= $(shell command -v include-what-you-use 2>/dev/null || echo include-what-you-use)
IWYU_TOOL      ?= $(shell command -v iwyu_tool 2>/dev/null || echo iwyu_tool)
FIX_INCLUDES   ?= $(shell command -v fix_includes.py 2>/dev/null || echo fix_includes.py)

# Default scope and options
IWYU_SCOPE       ?= $(SRCDIR)
IWYU_TOOL_ARGS   ?= -j $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) -p .
IWYU_IWYU_ARGS   ?= -Xiwyu --no_fwd_decls
IWYU_OUT         ?= iwyu.out

# Capture tool (bear preferred)
BEAR             := $(shell command -v bear 2>/dev/null)
INTERCEPT        := $(shell command -v intercept-build 2>/dev/null)

.PHONY: iwyu iwyu-apply iwyu-diff iwyu-file iwyu-clean compile_commands.json

compile_commands.json:
	@if [ -n "$(BEAR)" ]; then \
	    echo "[IWYU] capturing CPU build compile commands with bear"; \
	    $(BEAR) -- $(MAKE) -B cpu; \
	elif [ -n "$(INTERCEPT)" ]; then \
	    echo "[IWYU] capturing CPU build compile commands with intercept-build"; \
	    $(INTERCEPT) --cdb compile_commands.json -- $(MAKE) -B cpu; \
	else \
	    echo "ERROR: need 'bear' or 'intercept-build' to produce compile_commands.json"; \
	    echo "       install with: sudo apt install bear"; \
	    exit 2; \
	fi; \
	if [ ! -f compile_commands.json ]; then \
	    echo "compile_commands.json not generated."; \
	    exit 2; \
	fi

iwyu: compile_commands.json
	@echo "[IWYU] analyzing CPU build (MKL=$(mkl), OpenMP=$(openmp))"
	@if $(IWYU_TOOL) $(IWYU_TOOL_ARGS) $(IWYU_SCOPE) -- $(IWYU_IWYU_ARGS) > "$(IWYU_OUT)" 2>&1; then \
	    echo "[IWYU] wrote $(IWYU_OUT)"; \
	    echo "Tip:  make iwyu-diff   # preview edits"; \
	    echo "      make iwyu-apply  # apply edits"; \
	else \
	    echo "[IWYU] finished with non-zero exit (check $(IWYU_OUT))"; \
	fi

iwyu-diff:
	@if [ ! -f "$(IWYU_OUT)" ]; then \
	    echo "Run 'make iwyu' first."; \
	    exit 2; \
	fi; \
	$(FIX_INCLUDES) --nosafe_headers --reorder --nocomments --dry_run < "$(IWYU_OUT)" | head -n 200

iwyu-apply:
	@if [ ! -f "$(IWYU_OUT)" ]; then \
	    echo "Run 'make iwyu' first."; \
	    exit 2; \
	fi; \
	$(FIX_INCLUDES) --nosafe_headers --reorder --nocomments < "$(IWYU_OUT)"; \
	echo "[IWYU] includes updated. Rebuild to verify."

iwyu-file: compile_commands.json
	@if [ -z "$(FILE)" ]; then \
	    echo "Usage: make iwyu-file FILE=src/your_file.cpp"; \
	    exit 2; \
	fi; \
	echo "[IWYU] analyzing $(FILE)"; \
	if $(IWYU_TOOL) $(IWYU_TOOL_ARGS) $(FILE) -- $(IWYU_IWYU_ARGS) > "iwyu.$(notdir $(FILE)).out" 2>&1; then \
	    echo "[IWYU] wrote iwyu.$(notdir $(FILE)).out"; \
	    echo "Apply: $(FIX_INCLUDES) --nosafe_headers --reorder --nocomments < iwyu.$(notdir $(FILE)).out"; \
	else \
	    echo "[IWYU] finished with non-zero exit (check iwyu.$(notdir $(FILE)).out)"; \
	fi

iwyu-clean:
	@rm -f "$(IWYU_OUT)" compile_commands.json iwyu.*.out
	@echo "[IWYU] cleaned"


#===============================================================
# Clean Targets
#===============================================================

clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib

