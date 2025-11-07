# ===============================================================
# FEMaster Makefile (Linux/WSL/macOS) — CPU/GPU, MKL, static OK
# ===============================================================

PROJECT_NAME := FEMaster
CXX          := g++
NVCC         := nvcc

# ---------------------------------------------------------------
# Feature toggles (override via CLI, e.g. make mkl=1 static_link=1)
# ---------------------------------------------------------------
openmp             ?= 1   # 1: enable OpenMP (-fopenmp)
mkl                ?= 1   # 1: link MKL
mkl_sequential     ?= 0   # 1: libmkl_sequential; 0: libmkl_intel_thread (+iomp5)
debug              ?= 0   # 1: -g3 and O2; 0: -DNDEBUG -DEIGEN_NO_DEBUG and O3
double_precision   ?= 1   # 1: -DDOUBLE_PRECISION (and CUDA_DOUBLE_PRECISION)
eigen_fast_compile ?= 1   # 1: Eigen compile-speed tweaks
time_report        ?= 0   # 1: -ftime-report
static_link        ?= 0   # 1: -static-libstdc++ -static-libgcc + MKL static archives
cuda_dp            ?= 1   # keep for legacy guards in GPU code paths

# -------- Normalize toggles (empty env vars override ?=) --------
define _norm_toggle
$1 := $(strip $($1))
ifeq ($($1),)
  $1 := $2
endif
endef
$(eval $(call _norm_toggle,openmp,1))
$(eval $(call _norm_toggle,mkl,1))
$(eval $(call _norm_toggle,mkl_sequential,0))
$(eval $(call _norm_toggle,debug,0))
$(eval $(call _norm_toggle,double_precision,1))
$(eval $(call _norm_toggle,eigen_fast_compile,1))
$(eval $(call _norm_toggle,time_report,0))
$(eval $(call _norm_toggle,static_link,0))
$(eval $(call _norm_toggle,cuda_dp,1))

# if sequential is requested, force MKL on
mkl := $(if $(filter 1,$(mkl_sequential)),1,$(mkl))

# ---------------------------------------------------------------
# Paths
# ---------------------------------------------------------------
SRCDIR   := src
TESTDIR  := tests
OBJDIR   := obj
BINDIR   := bin
TARGET   := $(BINDIR)/$(PROJECT_NAME)

# Adjust /usr/local/include if Eigen/argparse live elsewhere
INCLUDES := -I/usr/local/include $(if $(MKLROOT),-I$(MKLROOT)/include)

# MKL paths
MKL_PATH := $(MKLROOT)/lib/intel64
ifeq ($(shell uname -s),Darwin)
  MKL_PATH := $(MKLROOT)/lib
endif

# Static Intel OpenMP (iomp5) archive (adjust if different on your system)
IOMP5_A  := /opt/intel/oneapi/compiler/2024.2/lib/libiomp5.a

# ---------------------------------------------------------------
# Sources / Objects
# ---------------------------------------------------------------
MAIN_SRC     := $(SRCDIR)/main.cpp
CPP_SRCS     := $(filter-out $(MAIN_SRC), $(sort $(shell find $(SRCDIR) -name '*.cpp')))
CU_SRCS      := $(sort $(shell find $(SRCDIR) -name '*.cu'))
TST_SRCS     := $(sort $(shell find $(TESTDIR) -name '*.cpp'))

GPP_OBJDIR   := $(OBJDIR)/gpp
NVCC_OBJROOT := $(OBJDIR)/nvcc
CPP_OBJDIR   := $(NVCC_OBJROOT)/cpp
CU_OBJDIR    := $(NVCC_OBJROOT)/cuda

GPP_OBJS     := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)
MAIN_OBJ_CPU := $(MAIN_SRC:$(SRCDIR)/%.cpp=$(GPP_OBJDIR)/%.o)

GPU_CPP_OBJS := $(CPP_SRCS:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
MAIN_OBJ_GPU := $(MAIN_SRC:$(SRCDIR)/%.cpp=$(CPP_OBJDIR)/%.o)
CU_OBJS      := $(CU_SRCS:$(SRCDIR)/%.cu=$(CU_OBJDIR)/%.o)

TST_OBJS     := $(TST_SRCS:$(TESTDIR)/%.cpp=$(GPP_OBJDIR)/%.o) $(GPP_OBJS)

# ---------------------------------------------------------------
# Tools
# ---------------------------------------------------------------
ifneq (,$(shell which ccache 2>/dev/null))
  CXX  := ccache $(CXX)
  NVCC := ccache $(NVCC)
endif

# ---------------------------------------------------------------
# Flags
# ---------------------------------------------------------------
CXXFLAGS := -std=c++17 $(INCLUDES) -MMD -MP
NVCCFLAGS := -std=c++17 $(INCLUDES) -MMD -MP --expt-relaxed-constexpr
FEATURES :=

# Debug vs Release
ifeq ($(debug),1)
  CXXFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(CXXFLAGS)) -O2 -g3
  NVCCFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(NVCCFLAGS)) -O2 -G -g
else
  CXXFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(CXXFLAGS)) -O3 -DNDEBUG -DEIGEN_NO_DEBUG
  NVCCFLAGS := $(filter-out -O0 -O1 -O2 -O3,$(NVCCFLAGS)) -O3
endif

# Optional time report
ifeq ($(time_report),1)
  CXXFLAGS += -ftime-report
endif

# OpenMP (compile-time)
ifeq ($(openmp),1)
  CXXFLAGS  += -fopenmp
  NVCCFLAGS += -Xcompiler=-fopenmp -diag-suppress 20014
  FEATURES  += -DEIGEN_DONT_PARALLELIZE
endif

# Double precision (CPU & CUDA guards unified)
ifeq ($(double_precision),1)
  FEATURES += -DDOUBLE_PRECISION -DCUDA_DOUBLE_PRECISION
endif

# Keep legacy cuda_dp guard if used in code
ifeq ($(cuda_dp),1)
  FEATURES += -DCUDA_DOUBLE_PRECISION
endif

# Eigen fast compile (keeps runtime same)
ifeq ($(eigen_fast_compile),1)
  FEATURES += -DEIGEN_UNROLLING_LIMIT=50 -DEIGEN_DONT_INLINE
endif

# ---------------------------------------------------------------
# MKL linking (Linux/macOS)
# ---------------------------------------------------------------
LDLIBS  :=
LDFLAGS :=

ifeq ($(static_link),1)
  LDFLAGS += -static-libstdc++ -static-libgcc
endif

ifeq ($(mkl),1)
  FEATURES += -DMKL_LP64 -DUSE_MKL
  # static vs dynamic MKL + runtime
  ifeq ($(static_link),1)
    ifeq ($(mkl_sequential),1)
      MKL_LIBS := -Wl,--start-group \
                    $(MKL_PATH)/libmkl_intel_lp64.a \
                    $(MKL_PATH)/libmkl_core.a \
                    $(MKL_PATH)/libmkl_sequential.a \
                  -Wl,--end-group \
                  -lpthread -lm -ldl
    else
      MKL_LIBS := -Wl,--start-group \
                    $(MKL_PATH)/libmkl_intel_lp64.a \
                    $(MKL_PATH)/libmkl_core.a \
                    $(MKL_PATH)/libmkl_intel_thread.a \
                  -Wl,--end-group \
                  -Wl,--whole-archive $(IOMP5_A) -Wl,--no-whole-archive \
                  -lpthread -lm -ldl
    endif
  else
    ifeq ($(mkl_sequential),1)
      MKL_LIBS := -L$(MKL_PATH) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
    else
      MKL_LIBS := -L$(MKL_PATH) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl
    endif
  endif
  LDLIBS += $(MKL_LIBS)
endif

# ---------------------------------------------------------------
# Link flags: avoid pulling libgomp when using iomp5
# ---------------------------------------------------------------
CXXFLAGS_LINK := $(CXXFLAGS)
ifeq ($(mkl),1)
  ifeq ($(mkl_sequential),0)
    CXXFLAGS_LINK := $(filter-out -fopenmp,$(CXXFLAGS))
  endif
endif

# ===============================================================
# Executables
# ===============================================================
EXE_CPU  := $(BINDIR)/$(PROJECT_NAME)
EXE_GPU  := $(BINDIR)/$(PROJECT_NAME)_gpu
EXE_TST  := $(BINDIR)/$(PROJECT_NAME)_test

# ===============================================================
# Targets
# ===============================================================
.PHONY: all cpu gpu tests clean info help pp-startup iwyu iwyu-apply iwyu-diff iwyu-file iwyu-clean compile_commands.json clean-exp-lib

all: info cpu gpu tests

cpu: info $(EXE_CPU)
gpu: CXXFLAGS  += -DSUPPORT_GPU
gpu: NVCCFLAGS += -DSUPPORT_GPU
gpu: info $(EXE_GPU) clean-exp-lib

tests: info $(EXE_TST)

# ---------------------------------------------------------------
# Link rules
# ---------------------------------------------------------------
$(EXE_CPU): $(GPP_OBJS) $(MAIN_OBJ_CPU)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS_LINK) $(FEATURES) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(EXE_GPU): $(GPU_CPP_OBJS) $(CU_OBJS) $(MAIN_OBJ_GPU)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(FEATURES) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(EXE_TST): $(TST_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(FEATURES) $^ $(LDFLAGS) -L/usr/local/lib -lgtest -lgtest_main -pthread -o $@

# ---------------------------------------------------------------
# Compile rules with timing
# ---------------------------------------------------------------
$(GPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(CXX) $(CXXFLAGS) $(FEATURES) -c $< -o $@; \
	END=$$(date +%s%N); \
	DUR=$$(echo "scale=3; ($$END-$$START)/1000000000" | bc); \
	printf "%-50s : %6s s\n" "Compiling $<" $$(echo $$DUR | awk '{printf "%.1f", $$1}')

$(CPP_OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(NVCC) $(NVCCFLAGS) $(FEATURES) -x cu -c $< -o $@; \
	END=$$(date +%s%N); \
	DUR=$$(echo "scale=3; ($$END-$$START)/1000000000" | bc); \
	printf "%-50s : %6s s\n" "Compiling $< (nvcc -x cu)" $$(echo $$DUR | awk '{printf "%.1f", $$1}')

$(CU_OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(NVCC) $(NVCCFLAGS) $(FEATURES) -c $< -o $@; \
	END=$$(date +%s%N); \
	DUR=$$(echo "scale=3; ($$END-$$START)/1000000000" | bc); \
	printf "%-50s : %6s s\n" "Compiling $< (cuda)" $$(echo $$DUR | awk '{printf "%.1f", $$1}')

$(GPP_OBJDIR)/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(@D)
	@START=$$(date +%s%N); \
	$(CXX) $(CXXFLAGS) $(FEATURES) -c $< -o $@; \
	END=$$(date +%s%N); \
	DUR=$$(echo "scale=3; ($$END-$$START)/1000000000" | bc); \
	printf "%-50s : %6s s\n" "Compiling test $<" $$(echo $$DUR | awk '{printf "%.1f", $$1}')

# deps
-include $(GPP_OBJS:.o=.d)
-include $(GPU_CPP_OBJS:.o=.d)
-include $(CU_OBJS:.o=.d)
-include $(MAIN_OBJ_CPU:.o=.d)
-include $(MAIN_OBJ_GPU:.o=.d)
-include $(TST_OBJS:.o=.d)

# ===============================================================
# Info / Utils
# ===============================================================
info:
	@echo "Build Configuration:"
	@echo "  Target           : $(TARGET)"
	@echo "  CXX              : $(CXX)"
	@echo "  NVCC             : $(NVCC)"
	@echo "  CXXFLAGS         : $(CXXFLAGS)"
	@echo "  NVCCFLAGS        : $(NVCCFLAGS)"
	@echo "  CXXFLAGS_LINK    : $(CXXFLAGS_LINK)"
	@echo "  FEATURES         : $(FEATURES)"
	@echo "  LDFLAGS          : $(LDFLAGS)"
	@echo "  LDLIBS           : $(LDLIBS)"
	@echo "  MKL Enabled      : $(mkl)"
	@echo "  MKL Sequential   : $(mkl_sequential)"
	@echo "  OpenMP           : $(openmp)"
	@echo "  Debug            : $(debug)"
 	@echo "  Double Precision : $(double_precision)"
	@echo "  Eigen Fast Comp. : $(eigen_fast_compile)"
	@echo "  Static Link      : $(static_link)"
	@echo "  MKLROOT          : $(MKLROOT)"
	@echo "  MKL_PATH         : $(MKL_PATH)"
	@echo "  IOMP5_A          : $(IOMP5_A)"

pp-startup:
	@echo "Preprocessor flags visible in startup.cpp (selected):"
	@$(CXX) $(CXXFLAGS) $(FEATURES) -dM -E $(SRCDIR)/core/startup.cpp | \
	  grep -E 'DOUBLE_PRECISION|CUDA_DOUBLE|USE_MKL|_OPENMP' || true

clean:
	@echo "Cleaning..."
	@rm -rf $(OBJDIR) $(BINDIR)

help:
	@echo "Usage examples:"
	@echo "  make -j                                 # default: MKL threaded, OpenMP on, release"
	@echo "  make -j mkl=1 mkl_sequential=1          # MKL sequential (keine iomp5-Abhängigkeit)"
	@echo "  make -j openmp=0                        # OpenMP komplett aus"
	@echo "  make -j debug=1                         # Debug (O2 + -g3)"
	@echo "  make -j static_link=1                   # statische libstdc++/libgcc + MKL-Archive"
	@echo "  make -j static_link=1 mkl_sequential=1  # komplett ohne iomp5 (sequentielles MKL)"
	@echo "  make clean"
	@echo "  make info"

# ===============================================================
# IWYU (Include-What-You-Use) — optional
# ===============================================================
IWYU           ?= $(shell command -v include-what-you-use 2>/dev/null || echo include-what-you-use)
IWYU_TOOL      ?= $(shell command -v iwyu_tool 2>/dev/null || echo iwyu_tool)
FIX_INCLUDES   ?= $(shell command -v fix_includes.py 2>/dev/null || echo fix_includes.py)

IWYU_SCOPE       ?= $(SRCDIR)
IWYU_TOOL_ARGS   ?= -j $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) -p .
IWYU_IWYU_ARGS   ?= -Xiwyu --no_fwd_decls
IWYU_OUT         ?= iwyu.out

BEAR             := $(shell command -v bear 2>/dev/null)
INTERCEPT        := $(shell command -v intercept-build 2>/dev/null)

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

# ===============================================================
# Clean extras
# ===============================================================
clean-exp-lib:
	@-rm -f $(BINDIR)/*.exp $(BINDIR)/*.lib
