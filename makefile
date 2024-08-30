# Compiler options
NVCC = nvcc
CXX = c++

CXXFLAGS = -std=c++17 -O3 -I ./include/
NVCCFLAGS = -std=c++17 -O3 -I ./include/ --expt-relaxed-constexpr
NVCCLIBS := -lcusolver -lcublas -lcusparse

UNAME := $(shell uname)
IS_CLANG := $(shell $(CXX) --version | grep -q "clang"; echo $$?)


# Directories
SRCDIR = src
TESTDIR = tests
OBJDIR = obj/
GPP_OBJDIR = $(OBJDIR)gpp/
NVCC_OBJDIR = $(OBJDIR)nvcc/
CPP_OBJDIR = $(NVCC_OBJDIR)cpp/
CU_OBJDIR  = $(NVCC_OBJDIR)/cuda/
BINDIR = bin

# Source files
MAIN_SRC := $(SRCDIR)/main.cpp
CPP_SRCS := $(sort $(filter-out $(MAIN_SRC), $(shell find $(SRCDIR) -name '*.cpp')))
CU_SRCS  := $(sort $(shell find $(SRCDIR) -name '*.cu'))
TST_SRCS := $(sort $(shell find $(TESTDIR) -name '*.cpp'))

# Object files for normal compilation
GPP_MAIN_OBJ := $(GPP_OBJDIR)/main.o
NVCC_MAIN_OBJ := $(CPP_OBJDIR)/main.o
GPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
GPU_CPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
CU_OBJS  := $(CU_SRCS:$(SRCDIR)/%.cu=$(CU_OBJDIR)/%.o)

# Test objects, excluding main.cpp
TST_OBJS := $(TST_SRCS:$(TESTDIR)/%.cpp=$(GPP_OBJDIR)/%.o) $(GPP_OBJS)

EXE_CPU  := $(BINDIR)/FEMaster.exe
EXE_GPU  := $(BINDIR)/FEMaster_gpu.exe
EXE_TST  := $(BINDIR)/test_suite.exe

# Optional flags
ifeq ($(cuda_double_precision), 1)
	CXXFLAGS  += -DCUDA_DOUBLE_PRECISION
	NVCCFLAGS += -DCUDA_DOUBLE_PRECISION
endif

ifeq ($(show_array_processes), 1)
	CXXFLAGS  += -DSHOW_ARRAY_PROCESSES
	NVCCFLAGS += -DSHOW_ARRAY_PROCESSES
endif

ifeq ($(debug), 1)
	CXXFLAGS  += -DNDEBUG
	NVCCFLAGS += -DNDEBUG
endif

all: cpu gpu tests

cpu: $(EXE_CPU)

gpu: CXXFLAGS += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: $(EXE_GPU) clean-exp-lib

tests: $(EXE_TST)

$(EXE_CPU): $(GPP_MAIN_OBJ) $(GPP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(EXE_GPU): $(NVCC_MAIN_OBJ) $(GPU_CPP_OBJS) $(CU_OBJS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) $^ -o $@

$(EXE_TST): $(TST_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -lgtest -lgtest_main -pthread -o $@

# Compilation rules for each type of object
$(GPP_OBJDIR)/main.o: $(SRCDIR)/main.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CPP_OBJDIR)/main.o: $(SRCDIR)/main.cpp
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -x cu -c $< -o $@

$(GPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -x cu -c $< -o $@

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(NVCCLIBS) -c $< -o $@

$(GPP_OBJDIR)/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib

clean:
	rm -rf $(GPP_OBJDIR) $(NVCC_OBJDIR) $(BINDIR) $(OBJDIR)
