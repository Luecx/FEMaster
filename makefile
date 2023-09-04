# Compiler options
NVCC = nvcc
CXX = g++
CXX = /opt/homebrew/opt/llvm/bin/clang++

CXXFLAGS = -std=c++17 -I ./include/ -fopenmp
NVCCFLAGS = -std=c++17 -O3 -I ./include/ --expt-relaxed-constexpr
NVCCLIBS := -lcusolver -lcublas -lcusparse

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CXXFLAGS += -fopenmp
	NVCCFLAGS += -Xcompiler=-fopenmp
endif

# Directories
SRCDIR = src
OBJDIR = obj/
GPP_OBJDIR = $(OBJDIR)gpp/
NVCC_OBJDIR = $(OBJDIR)nvcc/
CPP_OBJDIR = $(NVCC_OBJDIR)cpp/
CU_OBJDIR  = $(NVCC_OBJDIR)cuda/
BINDIR = bin

# Files
CPP_SRCS := $(sort $(shell find $(SRCDIR) -name '*.cpp'))
CU_SRCS  := $(sort $(shell find $(SRCDIR) -name '*.cu'))

GPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
GPU_CPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
CU_OBJS  := $(CU_SRCS:$(SRCDIR)/%.cu=$(CU_OBJDIR)/%.o)

EXE_CPU  := $(BINDIR)/FEMaster.exe
EXE_GPU  := $(BINDIR)/FEMaster_gpu.exe

# Check if double_precision is set to 1
ifeq ($(double_precision), 1)
	CXXFLAGS  += -DDOUBLE_PRECISION
	NVCCFLAGS += -DDOUBLE_PRECISION
endif

# check if show_array_processes is set to 1
ifeq ($(show_array_processes), 1)
	CXXFLAGS  += -DSHOW_ARRAY_PROCESSES
	NVCCFLAGS += -DSHOW_ARRAY_PROCESSES
endif

# check if debug is set to 1
ifeq ($(debug), 1)
	CXXFLAGS  += -DNDEBUG
	NVCCFLAGS += -DNDEBUG
endif

all: cpu gpu

cpu: $(EXE_CPU)

gpu: CXXFLAGS += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: $(EXE_GPU) clean-exp-lib

$(EXE_CPU): $(GPP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $(EXE_CPU)

$(EXE_GPU): $(GPU_CPP_OBJS) $(CU_OBJS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $(GPU_CPP_OBJS) $(CU_OBJS) -o $(EXE_GPU)

$(GPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -x cu -c $< -o $@

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) -c $< -o $@

clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib

clean:
	rm -rf $(GPP_OBJDIR) $(NVCC_OBJDIR) $(BINDIR) $(OBJDIR)
