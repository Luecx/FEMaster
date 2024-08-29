#===============================================================
# Compiler Options and Configurations
#===============================================================

# Compiler settings
NVCC = nvcc
CXX  = g++

#===============================================================
# Default Feature Flags
#===============================================================

openmp   ?= 1   # Enable/disable OpenMP
mkl      ?= 0   # Enable/disable MKL (Math Kernel Library)
cuda_dp  ?= 1   # Enable/disable CUDA Double Precision
ar_pcs   ?= 0   # Enable/disable Show Array Processes
debug    ?= 0   # Enable/disable Debug mode

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

# MKL support
ifeq ($(mkl), 1)
    FEATURE_FLAGS += -DMKL_LP64
    FEATURE_FLAGS += -DUSE_MKL
    LIBS += $(MKL_LIBS)
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
endif

#===============================================================
# MKL Library Paths and Libraries
#===============================================================

MKLROOT   ?= $(shell echo ${MKLROOT})

MKL_LIBS := \
    -Wl,--start-group \
    $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
    $(MKLROOT)/lib/intel64/libmkl_sequential.a \
    $(MKLROOT)/lib/intel64/libmkl_core.a \
    -Wl,--end-group \
    -lpthread \
    -lm \
    -ldl

#===============================================================
# CUDA Libraries
#===============================================================

NVCCLIBS := \
    -lcusolver \
    -lcublas \
    -lcusparse

#===============================================================
# System Information (Optional Flags Based on System)
#===============================================================

UNAME := $(shell uname)

#===============================================================
# Directories and File Paths
#===============================================================

SRCDIR      = src
OBJDIR      = obj/
GPP_OBJDIR  = $(OBJDIR)gpp/
NVCC_OBJDIR = $(OBJDIR)nvcc/
CPP_OBJDIR  = $(NVCC_OBJDIR)cpp/
CU_OBJDIR   = $(NVCC_OBJDIR)cuda/
BINDIR      = bin

#===============================================================
# Source and Object Files
#===============================================================

CPP_SRCS := $(sort $(shell find $(SRCDIR) -name '*.cpp'))
CU_SRCS  := $(sort $(shell find $(SRCDIR) -name '*.cu'))

GPP_OBJS     := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
GPU_CPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
CU_OBJS      := $(CU_SRCS:$(SRCDIR)/%.cu=$(CU_OBJDIR)/%.o)

#===============================================================
# Executable Files
#===============================================================

EXE_CPU  := $(BINDIR)/FEMaster.exe
EXE_GPU  := $(BINDIR)/FEMaster_gpu.exe

#===============================================================
# Display Compilation Flags
#===============================================================
#===============================================================
# Display Compilation Flags
#===============================================================

show_flags:
	@echo "Compilation Flags:"
	@echo "-------------------"
	@echo "C++ Compiler : $(CXX)"
	@echo "NVCC Compiler: $(NVCC)"
	@echo "CXXFLAGS     : $(CXXFLAGS)"
	@echo "NVCCFLAGS    : $(NVCCFLAGS)"
	@echo "LIBS         :"
	@echo "$(LIBS)" | tr ' ' '\n' | sed 's/^/    /'
	@echo "Feature Flags: $(FEATURE_FLAGS)"
	@echo ""
	@echo "Feature Flags:"
	@echo "-------------------"
	@echo "OpenMP               : $(if $(findstring -fopenmp, $(CXXFLAGS)),Enabled,Disabled)"
	@echo "MKL                  : $(if $(findstring -DUSE_MKL, $(FEATURE_FLAGS)),Enabled,Disabled)"
	@echo "CUDA Double Precision: $(if $(findstring -DCUDA_DOUBLE_PRECISION, $(FEATURE_FLAGS)),Enabled,Disabled)"
	@echo "Show Array Processes : $(if $(findstring -DSHOW_ARRAY_PROCESSES, $(FEATURE_FLAGS)),Enabled,Disabled)"
	@echo "Debug                : $(if $(findstring -DNDEBUG, $(FEATURE_FLAGS)),Enabled,Disabled)"


#===============================================================
# Build Targets
#===============================================================

all: show_flags cpu gpu

cpu: show_flags $(EXE_CPU)

gpu: CXXFLAGS  += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: show_flags $(EXE_GPU) clean-exp-lib

#===============================================================
# Build Rules
#===============================================================

CXXFLAGS   += -MMD -MP
NVCCFLAGS  += -MMD -MP

$(EXE_CPU): $(GPP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) $^ $(LIBS) -o $@

$(EXE_GPU): $(GPU_CPP_OBJS) $(CU_OBJS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $(GPU_CPP_OBJS) $(CU_OBJS) $(LIBS) -o $@

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

