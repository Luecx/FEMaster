# Experimental AMGX integration for the FEG/amgx benchmark branch.
# Load with:
#   cmake ... -DCMAKE_PROJECT_INCLUDE="$PWD/cmake/amgx.cmake"

set(FEMASTER_AMGX_ROOT "$ENV{AMGX_ROOT}" CACHE PATH
        "Root of a local NVIDIA AMGX source/build tree"
)

if(NOT FEMASTER_ENABLE_CUDA)
    message(FATAL_ERROR "AMGX requires FEMASTER_ENABLE_CUDA=ON")
endif()

set(FEMASTER_AMGX_HINTS)
if(FEMASTER_AMGX_ROOT)
    list(APPEND FEMASTER_AMGX_HINTS "${FEMASTER_AMGX_ROOT}")
endif()
if(DEFINED ENV{AMGX_ROOT})
    list(APPEND FEMASTER_AMGX_HINTS "$ENV{AMGX_ROOT}")
endif()

find_path(FEMASTER_AMGX_INCLUDE_DIR
        NAMES amgx_c.h
        HINTS ${FEMASTER_AMGX_HINTS}
        PATH_SUFFIXES include
)

find_library(FEMASTER_AMGX_LIBRARY
        NAMES amgxsh
        HINTS ${FEMASTER_AMGX_HINTS}
        PATH_SUFFIXES build lib lib64
)

if(NOT FEMASTER_AMGX_INCLUDE_DIR)
    message(FATAL_ERROR
            "AMGX header amgx_c.h was not found. Set AMGX_ROOT or FEMASTER_AMGX_ROOT."
    )
endif()

if(NOT FEMASTER_AMGX_LIBRARY)
    message(FATAL_ERROR
            "AMGX shared library libamgxsh.so was not found. Build AMGX first and set AMGX_ROOT or FEMASTER_AMGX_ROOT."
    )
endif()

# This file is injected at the end of project(), before FEMaster creates targets,
# so directory-level settings apply to femaster_core and FEMaster.
include_directories(SYSTEM "${FEMASTER_AMGX_INCLUDE_DIR}")
link_libraries("${FEMASTER_AMGX_LIBRARY}")
add_compile_definitions(USE_AMGX)

get_filename_component(FEMASTER_AMGX_LIBRARY_DIR "${FEMASTER_AMGX_LIBRARY}" DIRECTORY)
list(PREPEND CMAKE_BUILD_RPATH "${FEMASTER_AMGX_LIBRARY_DIR}")

message(STATUS "AMGX include : ${FEMASTER_AMGX_INCLUDE_DIR}")
message(STATUS "AMGX library : ${FEMASTER_AMGX_LIBRARY}")
