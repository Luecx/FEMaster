cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0091 NEW)

project(FEMaster VERSION 2.0.0 LANGUAGES C CXX)

include(FetchContent)

# -----------------------------------------------------------------------------
# User options
# -----------------------------------------------------------------------------
option(FEMASTER_BUILD_TESTS "Build and register the GoogleTest executable" OFF)
option(FEMASTER_ENABLE_CUDA "Build CUDA solver backends" OFF)
option(FEMASTER_ENABLE_CUDSS "Use NVIDIA cuDSS for CUDA sparse direct solves" OFF)
option(FEMASTER_ENABLE_MKL "Use Intel oneMKL for Eigen operations and PARDISO" OFF)
option(FEMASTER_ENABLE_ACCELERATE "Use Apple Accelerate for sparse direct solves" ${APPLE})
option(FEMASTER_MKL_SEQUENTIAL "Use the sequential oneMKL threading layer" OFF)
option(FEMASTER_MKL_STATIC "Link the oneMKL component libraries statically" ON)
option(FEMASTER_ENABLE_OPENMP "Enable OpenMP for FEMaster CPU loops" OFF)
option(FEMASTER_DOUBLE_PRECISION "Use double precision arithmetic" ON)
option(FEMASTER_FAST_EIGEN_COMPILE "Use Eigen compile-time reduction flags" ON)
option(FEMASTER_STATIC_RUNTIME "Link compiler C/C++ runtimes statically where supported" OFF)
option(FEMASTER_FULLY_STATIC "Produce a fully static Linux executable" OFF)
option(FEMASTER_TIME_REPORT "Enable compiler timing diagnostics where supported" OFF)

set(FEMASTER_CUDSS_ROOT "$ENV{CUDSS_DIR}" CACHE PATH "Root of an NVIDIA cuDSS installation")

set(FEMASTER_EIGEN_VERSION "5.0.1")
set(FEMASTER_SPECTRA_VERSION "v1.2.0")
set(FEMASTER_ARGPARSE_VERSION "v3.2")
set(FEMASTER_LZ4_VERSION "v1.10.0")
set(FEMASTER_GTEST_VERSION "v1.18.0")

set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/_deps" CACHE PATH "FetchContent dependency directory")

# -----------------------------------------------------------------------------
# Global build behavior
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)

# Select the runtime before creating targets so LZ4 and GoogleTest match FEMaster.
if(MSVC AND FEMASTER_STATIC_RUNTIME)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

if(CMAKE_CONFIGURATION_TYPES)
    foreach(config Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${config}" config_upper)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin/${config}")
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/lib/${config}")
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/lib/${config}")
    endforeach()
endif()

# -----------------------------------------------------------------------------
# Option validation
# -----------------------------------------------------------------------------
if(FEMASTER_ENABLE_CUDSS AND NOT FEMASTER_ENABLE_CUDA)
    message(FATAL_ERROR "FEMASTER_ENABLE_CUDSS requires FEMASTER_ENABLE_CUDA=ON")
endif()

if(FEMASTER_ENABLE_CUDA AND APPLE)
    message(FATAL_ERROR "CUDA is not supported by current NVIDIA toolkits on macOS")
endif()

if(FEMASTER_ENABLE_ACCELERATE AND NOT APPLE)
    message(FATAL_ERROR "FEMASTER_ENABLE_ACCELERATE is supported only on Apple platforms")
endif()

if(FEMASTER_ENABLE_ACCELERATE AND FEMASTER_ENABLE_MKL)
    message(FATAL_ERROR "FEMASTER_ENABLE_ACCELERATE and FEMASTER_ENABLE_MKL are mutually exclusive")
endif()

if(FEMASTER_FULLY_STATIC AND NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "FEMASTER_FULLY_STATIC is supported only on Linux")
endif()

if(FEMASTER_FULLY_STATIC AND FEMASTER_ENABLE_CUDA)
    message(FATAL_ERROR "FEMASTER_FULLY_STATIC is incompatible with CUDA shared libraries")
endif()

if(FEMASTER_FULLY_STATIC AND FEMASTER_ENABLE_MKL AND NOT FEMASTER_MKL_STATIC)
    message(FATAL_ERROR "FEMASTER_FULLY_STATIC requires FEMASTER_MKL_STATIC=ON")
endif()

if(FEMASTER_ENABLE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
endif()

# =============================================================================
# Dependencies
# =============================================================================

# Header-only dependencies. SOURCE_SUBDIR deliberately points to a nonexistent
# directory so their own CMakeLists.txt files are not executed.
FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG ${FEMASTER_EIGEN_VERSION}
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _femaster_no_cmake
)

FetchContent_Declare(spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG ${FEMASTER_SPECTRA_VERSION}
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _femaster_no_cmake
)

FetchContent_Declare(argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
        GIT_TAG ${FEMASTER_ARGPARSE_VERSION}
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _femaster_no_cmake
)

FetchContent_Declare(lz4
        GIT_REPOSITORY https://github.com/lz4/lz4.git
        GIT_TAG ${FEMASTER_LZ4_VERSION}
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR _femaster_no_cmake
)

FetchContent_MakeAvailable(eigen spectra argparse lz4)

add_library(femaster_eigen INTERFACE)
target_include_directories(femaster_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")
add_library(Eigen3::Eigen ALIAS femaster_eigen)

add_library(femaster_spectra INTERFACE)
target_include_directories(femaster_spectra SYSTEM INTERFACE "${spectra_SOURCE_DIR}/include")
target_link_libraries(femaster_spectra INTERFACE Eigen3::Eigen)
add_library(Spectra::Spectra ALIAS femaster_spectra)

add_library(femaster_argparse INTERFACE)
target_include_directories(femaster_argparse SYSTEM INTERFACE "${argparse_SOURCE_DIR}/include")
target_compile_features(femaster_argparse INTERFACE cxx_std_17)
add_library(argparse::argparse ALIAS femaster_argparse)

file(GLOB FEMASTER_LZ4_SOURCES CONFIGURE_DEPENDS "${lz4_SOURCE_DIR}/lib/*.c")
add_library(femaster_lz4 STATIC ${FEMASTER_LZ4_SOURCES})
target_include_directories(femaster_lz4 SYSTEM PUBLIC "${lz4_SOURCE_DIR}/lib")
target_compile_definitions(femaster_lz4 PRIVATE XXH_NAMESPACE=LZ4_)
set_target_properties(femaster_lz4 PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(LZ4::lz4 ALIAS femaster_lz4)

# =============================================================================
# FEMaster
# =============================================================================

file(GLOB_RECURSE FEMASTER_CPP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM FEMASTER_CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

set(FEMASTER_SOURCES ${FEMASTER_CPP_SOURCES})

if(FEMASTER_ENABLE_CUDA)
    file(GLOB FEMASTER_CUDA_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/cuda/*.cu")
    list(APPEND FEMASTER_SOURCES ${FEMASTER_CUDA_SOURCES})
endif()

add_library(femaster_core STATIC ${FEMASTER_SOURCES})
target_include_directories(femaster_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(femaster_core PUBLIC
        Eigen3::Eigen
        Spectra::Spectra
        argparse::argparse
        LZ4::lz4
)

target_compile_definitions(femaster_core PUBLIC
        $<$<BOOL:${FEMASTER_DOUBLE_PRECISION}>:DOUBLE_PRECISION>
        $<$<BOOL:${FEMASTER_DOUBLE_PRECISION}>:CUDA_DOUBLE_PRECISION>
        $<$<BOOL:${FEMASTER_FAST_EIGEN_COMPILE}>:EIGEN_UNROLLING_LIMIT=50>
        $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>:EIGEN_NO_DEBUG>
)

# -----------------------------------------------------------------------------
# Compiler options
# -----------------------------------------------------------------------------
set(FEMASTER_MSVC_MP_JOBS "8" CACHE STRING "MSVC /MP parallel compiler processes")

if(MSVC)
    target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/MP${FEMASTER_MSVC_MP_JOBS}>
            $<$<COMPILE_LANGUAGE:CXX>:/bigobj>
            $<$<COMPILE_LANGUAGE:CXX>:/permissive->
    )
else()
    target_compile_options(femaster_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>)
endif()

if(FEMASTER_TIME_REPORT)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(femaster_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ftime-report>)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(femaster_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ftime-trace>)
    elseif(MSVC)
        target_compile_options(femaster_core PRIVATE /Bt+)
    endif()
endif()

# -----------------------------------------------------------------------------
# Intel oneMKL
# -----------------------------------------------------------------------------
set(FEMASTER_MKL_LINK "disabled")
set(FEMASTER_MKL_THREADING "disabled")
set(FEMASTER_OPENMP_LINK_OPTIONS)

unset(FEMASTER_IOMP5_DLL CACHE)
unset(FEMASTER_IOMP5_DLL)

if(FEMASTER_ENABLE_MKL)
    set(MKL_INTERFACE "lp64")
    set(MKL_LINK "dynamic")
    set(MKL_THREADING "intel_thread")

    if(FEMASTER_MKL_STATIC)
        set(MKL_LINK "static")
    endif()

    if(FEMASTER_MKL_SEQUENTIAL)
        set(MKL_THREADING "sequential")
    endif()

    # For static threaded MKL, also select Intel's static OpenMP runtime.
    if(FEMASTER_MKL_STATIC AND NOT FEMASTER_MKL_SEQUENTIAL AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        find_file(FEMASTER_IOMP5_STATIC
                NAMES libiomp5.a
                HINTS
                "$ENV{CMPLR_ROOT}"
                "$ENV{ONEAPI_ROOT}/compiler/latest"
                "$ENV{MKLROOT}/../../compiler/latest"
                PATH_SUFFIXES lib compiler/lib linux/compiler/lib/intel64_lin
                DOC "Static Intel OpenMP runtime for threaded static oneMKL"
        )

        if(NOT FEMASTER_IOMP5_STATIC)
            message(FATAL_ERROR
                    "Static threaded oneMKL requires libiomp5.a. Set FEMASTER_IOMP5_STATIC, "
                    "enable FEMASTER_MKL_SEQUENTIAL, or disable FEMASTER_MKL_STATIC."
            )
        endif()

        set(OMP_LIBRARY "${FEMASTER_IOMP5_STATIC}")
    else()
        unset(OMP_LIBRARY CACHE)
        unset(OMP_LIBRARY)
    endif()

    find_package(MKL CONFIG REQUIRED
            HINTS "$ENV{MKLROOT}"
            PATH_SUFFIXES lib/cmake/mkl
    )

    target_link_libraries(femaster_core PUBLIC MKL::MKL)
    target_compile_definitions(femaster_core PUBLIC USE_MKL MKL_LP64)

    set(FEMASTER_MKL_LINK "${MKL_LINK}")
    set(FEMASTER_MKL_THREADING "${MKL_THREADING}")

    if(WIN32 AND NOT FEMASTER_MKL_SEQUENTIAL)
        find_file(FEMASTER_IOMP5_DLL
                NAMES libiomp5md.dll
                HINTS "${OMP_DLL_DIR}" "$ENV{ONEAPI_ROOT}/compiler/latest"
                PATH_SUFFIXES bin windows/redist/intel64_win/compiler
                DOC "Intel OpenMP runtime DLL copied beside FEMaster"
        )

        if(NOT FEMASTER_IOMP5_DLL)
            message(FATAL_ERROR "Threaded oneMKL requires libiomp5md.dll on Windows")
        endif()
    endif()
endif()

# -----------------------------------------------------------------------------
# Apple Accelerate
# -----------------------------------------------------------------------------
if(FEMASTER_ENABLE_ACCELERATE)
    find_library(FEMASTER_ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    target_link_libraries(femaster_core PUBLIC "${FEMASTER_ACCELERATE_FRAMEWORK}")
    target_compile_definitions(femaster_core PUBLIC USE_ACCELERATE)
endif()

# -----------------------------------------------------------------------------
# OpenMP
# -----------------------------------------------------------------------------
set(FEMASTER_OPENMP_RUNTIME "disabled")

if(FEMASTER_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_compile_definitions(femaster_core PUBLIC EIGEN_DONT_PARALLELIZE)

    message(STATUS "OpenMP_CXX_FOUND    : ${OpenMP_CXX_FOUND}")
    message(STATUS "OpenMP_CXX_FLAGS    : ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_CXX_LIB_NAMES: ${OpenMP_CXX_LIB_NAMES}")

    if(FEMASTER_ENABLE_MKL AND NOT FEMASTER_MKL_SEQUENTIAL)
        if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            message(FATAL_ERROR
                    "MSVC OpenMP and threaded oneMKL use different runtimes. "
                    "Disable OpenMP, select sequential MKL, or use IntelLLVM."
            )
        endif()

        if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
            target_compile_options(femaster_core PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/Qiopenmp>)
            list(APPEND FEMASTER_OPENMP_LINK_OPTIONS /Qiopenmp)
            set(FEMASTER_OPENMP_RUNTIME "Intel OpenMP via oneMKL (/Qiopenmp)")

        elseif(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(femaster_core PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/openmp:llvm>)
            list(APPEND FEMASTER_OPENMP_LINK_OPTIONS /openmp:llvm)
            set(FEMASTER_OPENMP_RUNTIME "LLVM OpenMP compile flag; runtime via oneMKL")

        else()
            if(OpenMP_CXX_FLAGS)
                separate_arguments(FEMASTER_OPENMP_FLAG_LIST NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
                target_compile_options(femaster_core PUBLIC ${FEMASTER_OPENMP_FLAG_LIST})
            endif()

            get_target_property(FEMASTER_OPENMP_COMPILE_OPTIONS
                    OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS
            )

            if(FEMASTER_OPENMP_COMPILE_OPTIONS)
                target_compile_options(femaster_core PUBLIC ${FEMASTER_OPENMP_COMPILE_OPTIONS})
            endif()

            set(FEMASTER_OPENMP_RUNTIME "Intel OpenMP via oneMKL")
        endif()
    else()
        target_link_libraries(femaster_core PUBLIC OpenMP::OpenMP_CXX)
        set(FEMASTER_OPENMP_RUNTIME "${OpenMP_CXX_LIB_NAMES}")
    endif()
endif()

# -----------------------------------------------------------------------------
# CUDA and cuDSS
# -----------------------------------------------------------------------------
set(FEMASTER_CUDA_VERSION "disabled")
set(FEMASTER_CUDSS_VERSION "disabled")

if(FEMASTER_ENABLE_CUDA)
    target_compile_definitions(femaster_core PUBLIC SUPPORT_GPU)
    target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    )

    target_link_libraries(femaster_core PUBLIC
            CUDA::cudart
            CUDA::cublas
            CUDA::cusparse
            CUDA::cusolver
    )

    set_target_properties(femaster_core PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED YES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    set(FEMASTER_CUDA_VERSION "${CUDAToolkit_VERSION}")

    if(FEMASTER_ENABLE_CUDSS)
        find_package(cudss CONFIG REQUIRED
                HINTS "${FEMASTER_CUDSS_ROOT}" "$ENV{CUDSS_DIR}"
                PATH_SUFFIXES lib/cmake/cudss
        )

        target_link_libraries(femaster_core PUBLIC cudss)
        target_compile_definitions(femaster_core PUBLIC USE_CUDSS)
        set(FEMASTER_CUDSS_VERSION "${cudss_VERSION}")
    endif()
endif()

# =============================================================================
# Executable
# =============================================================================
add_executable(FEMaster src/main.cpp)
target_link_libraries(FEMaster PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,femaster_core>")

if(FEMASTER_OPENMP_LINK_OPTIONS)
    target_link_options(FEMaster PRIVATE ${FEMASTER_OPENMP_LINK_OPTIONS})
endif()

if(FEMASTER_STATIC_RUNTIME AND NOT MSVC AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_options(FEMaster PRIVATE -static-libgcc -static-libstdc++)
endif()

if(FEMASTER_FULLY_STATIC)
    target_link_options(FEMaster PRIVATE -static)
endif()

# -----------------------------------------------------------------------------
# Runtime files
# -----------------------------------------------------------------------------
add_custom_command(TARGET FEMaster POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt"
        "$<TARGET_FILE_DIR:FEMaster>/LICENSE.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/THIRD_PARTY_NOTICES.txt"
        "$<TARGET_FILE_DIR:FEMaster>/THIRD_PARTY_NOTICES.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/licenses"
        "$<TARGET_FILE_DIR:FEMaster>/licenses"
        COMMENT "Copying license and third-party notice files"
)

if(FEMASTER_IOMP5_DLL)
    add_custom_command(TARGET FEMaster POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${FEMASTER_IOMP5_DLL}" "$<TARGET_FILE_DIR:FEMaster>"
            COMMENT "Copying Intel OpenMP runtime"
    )
endif()

if(WIN32 AND FEMASTER_ENABLE_CUDSS)
    add_custom_command(TARGET FEMaster POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:cudss>" "$<TARGET_FILE_DIR:FEMaster>"
            COMMENT "Copying cuDSS runtime"
    )
endif()

# =============================================================================
# Tests
# =============================================================================
if(FEMASTER_BUILD_TESTS)
    enable_testing()

    FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG ${FEMASTER_GTEST_VERSION}
            GIT_SHALLOW TRUE
    )

    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

    if(MSVC)
        if(FEMASTER_STATIC_RUNTIME)
            set(gtest_force_shared_crt OFF CACHE BOOL "" FORCE)
        else()
            set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        endif()
    endif()

    FetchContent_MakeAvailable(googletest)

    file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

    add_executable(FEMaster_test ${TEST_SOURCES})
    target_link_libraries(FEMaster_test PRIVATE
            "$<LINK_LIBRARY:WHOLE_ARCHIVE,femaster_core>"
            GTest::gtest
            GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(FEMaster_test)
endif()

# =============================================================================
# Configuration summary
# =============================================================================
if(CMAKE_CONFIGURATION_TYPES)
    set(FEMASTER_BUILD_MODE "multi-config: ${CMAKE_CONFIGURATION_TYPES}")
else()
    set(FEMASTER_BUILD_MODE "${CMAKE_BUILD_TYPE}")
endif()

message(STATUS "")
message(STATUS "================ FEMaster configuration ================")
message(STATUS "  Host system          : ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "  Target system        : ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  Build mode           : ${FEMASTER_BUILD_MODE}")
message(STATUS "  C++ compiler         : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Output directory     : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "  Dependencies         : ${FETCHCONTENT_BASE_DIR}")
message(STATUS "  Eigen                : ${FEMASTER_EIGEN_VERSION}")
message(STATUS "  Spectra              : ${FEMASTER_SPECTRA_VERSION}")
message(STATUS "  argparse             : ${FEMASTER_ARGPARSE_VERSION}")
message(STATUS "  LZ4                  : ${FEMASTER_LZ4_VERSION}")
message(STATUS "  Double precision     : ${FEMASTER_DOUBLE_PRECISION}")
message(STATUS "  OpenMP loops         : ${FEMASTER_ENABLE_OPENMP} (${FEMASTER_OPENMP_RUNTIME})")
message(STATUS "  oneMKL               : ${FEMASTER_ENABLE_MKL}")
message(STATUS "    link/threading     : ${FEMASTER_MKL_LINK} / ${FEMASTER_MKL_THREADING}")
message(STATUS "  Apple Accelerate     : ${FEMASTER_ENABLE_ACCELERATE}")
message(STATUS "  CUDA                 : ${FEMASTER_ENABLE_CUDA} (${FEMASTER_CUDA_VERSION})")
message(STATUS "  cuDSS                : ${FEMASTER_ENABLE_CUDSS} (${FEMASTER_CUDSS_VERSION})")
message(STATUS "  Static C/C++ runtime : ${FEMASTER_STATIC_RUNTIME}")
message(STATUS "  Fully static Linux   : ${FEMASTER_FULLY_STATIC}")
message(STATUS "  Tests                : ${FEMASTER_BUILD_TESTS}")
message(STATUS "==========================================================")
message(STATUS "")

add_custom_target(femaster-info
        COMMAND ${CMAKE_COMMAND} -E echo "FEMaster ${PROJECT_VERSION}"
        COMMAND ${CMAKE_COMMAND} -E echo "Target: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}"
        COMMAND ${CMAKE_COMMAND} -E echo "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
        COMMAND ${CMAKE_COMMAND} -E echo "Dependencies: ${FETCHCONTENT_BASE_DIR}"
        COMMAND ${CMAKE_COMMAND} -E echo "Eigen: ${FEMASTER_EIGEN_VERSION}"
        COMMAND ${CMAKE_COMMAND} -E echo "Spectra: ${FEMASTER_SPECTRA_VERSION}"
        COMMAND ${CMAKE_COMMAND} -E echo "OpenMP: ${FEMASTER_ENABLE_OPENMP} (${FEMASTER_OPENMP_RUNTIME})"
        COMMAND ${CMAKE_COMMAND} -E echo "oneMKL: ${FEMASTER_ENABLE_MKL} (${FEMASTER_MKL_LINK}, ${FEMASTER_MKL_THREADING})"
        COMMAND ${CMAKE_COMMAND} -E echo "Apple Accelerate: ${FEMASTER_ENABLE_ACCELERATE}"
        COMMAND ${CMAKE_COMMAND} -E echo "CUDA: ${FEMASTER_ENABLE_CUDA} (${FEMASTER_CUDA_VERSION})"
        COMMAND ${CMAKE_COMMAND} -E echo "cuDSS: ${FEMASTER_ENABLE_CUDSS} (${FEMASTER_CUDSS_VERSION})"
        VERBATIM
)
