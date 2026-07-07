cmake_minimum_required(VERSION 3.24)

cmake_policy(SET CMP0091 NEW)

project(FEMaster VERSION 2.0.0 LANGUAGES CXX)

# -----------------------------------------------------------------------------
# User options
# -----------------------------------------------------------------------------
option(FEMASTER_BUILD_TESTS "Build and register the GoogleTest executable" OFF)
option(FEMASTER_ENABLE_CUDA "Build CUDA solver backends" OFF)
option(FEMASTER_ENABLE_CUDSS "Use NVIDIA cuDSS for CUDA sparse direct solves" OFF)
option(FEMASTER_ENABLE_MKL "Use Intel oneMKL for Eigen operations and PARDISO" OFF)
option(FEMASTER_MKL_SEQUENTIAL "Use the sequential oneMKL threading layer" OFF)
option(FEMASTER_MKL_STATIC "Link the oneMKL component libraries statically" ON)
option(FEMASTER_ENABLE_OPENMP "Enable OpenMP for FEMaster CPU loops" OFF)
option(FEMASTER_DOUBLE_PRECISION "Use double precision arithmetic" ON)
option(FEMASTER_FAST_EIGEN_COMPILE "Use Eigen compile-time reduction flags" ON)
option(FEMASTER_STATIC_RUNTIME "Link compiler C/C++ runtimes statically where supported" OFF)
option(FEMASTER_FULLY_STATIC "Produce a fully static Linux executable" OFF)
option(FEMASTER_TIME_REPORT "Enable compiler timing diagnostics where supported" OFF)

set(FEMASTER_DEPS_INCLUDE_DIR "" CACHE PATH
    "Include root containing Eigen/, Spectra/, and argparse/"
)
set(FEMASTER_CUDSS_ROOT "$ENV{CUDSS_DIR}" CACHE PATH
    "Root of an NVIDIA cuDSS installation"
)

# -----------------------------------------------------------------------------
# Global build behavior
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel
    )
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

if(FEMASTER_ENABLE_CUDSS AND NOT FEMASTER_ENABLE_CUDA)
    message(FATAL_ERROR "FEMASTER_ENABLE_CUDSS requires FEMASTER_ENABLE_CUDA=ON")
endif()

if(FEMASTER_ENABLE_CUDA AND APPLE)
    message(FATAL_ERROR "CUDA is not supported by current NVIDIA toolkits on macOS")
endif()

if(FEMASTER_FULLY_STATIC AND NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux"))
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

# -----------------------------------------------------------------------------
# Header-only dependencies
# -----------------------------------------------------------------------------
set(FEMASTER_DEPS_INCLUDE_HINTS)
if(FEMASTER_DEPS_INCLUDE_DIR)
    list(APPEND FEMASTER_DEPS_INCLUDE_HINTS "${FEMASTER_DEPS_INCLUDE_DIR}")
endif()
if(DEFINED ENV{FEMASTER_DEPS_INCLUDE_DIR})
    list(APPEND FEMASTER_DEPS_INCLUDE_HINTS "$ENV{FEMASTER_DEPS_INCLUDE_DIR}")
endif()

if(WIN32)
    list(APPEND FEMASTER_DEPS_INCLUDE_HINTS
        "//wsl.localhost/Ubuntu-24.04/usr/local/include"
        "//wsl$/Ubuntu-24.04/usr/local/include"
    )
elseif(APPLE)
    list(APPEND FEMASTER_DEPS_INCLUDE_HINTS
        "/opt/homebrew/include"
        "/usr/local/include"
        "/usr/include"
    )
else()
    list(APPEND FEMASTER_DEPS_INCLUDE_HINTS
        "/usr/local/include"
        "/usr/include"
    )
endif()

find_package(Eigen3 CONFIG QUIET)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Core
        HINTS ${FEMASTER_DEPS_INCLUDE_HINTS}
        PATH_SUFFIXES eigen3
        DOC "Directory containing Eigen/"
    )
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Eigen was not found. Install its CMake package or set "
            "FEMASTER_DEPS_INCLUDE_DIR/EIGEN3_INCLUDE_DIR."
        )
    endif()
endif()

find_package(Spectra CONFIG QUIET)
if(NOT Spectra_FOUND)
    find_path(SPECTRA_INCLUDE_DIR
        NAMES Spectra/SymEigsSolver.h
        HINTS ${FEMASTER_DEPS_INCLUDE_HINTS}
        DOC "Directory containing Spectra/"
    )
    if(NOT SPECTRA_INCLUDE_DIR)
        message(FATAL_ERROR
            "Spectra was not found. Install its CMake package or set "
            "FEMASTER_DEPS_INCLUDE_DIR/SPECTRA_INCLUDE_DIR."
        )
    endif()
endif()

find_package(argparse CONFIG QUIET)
if(NOT argparse_FOUND)
    find_path(ARGPARSE_INCLUDE_DIR
        NAMES argparse/argparse.hpp
        HINTS ${FEMASTER_DEPS_INCLUDE_HINTS}
        DOC "Directory containing argparse/"
    )
    if(NOT ARGPARSE_INCLUDE_DIR)
        message(FATAL_ERROR
            "argparse was not found. Install its CMake package or set "
            "FEMASTER_DEPS_INCLUDE_DIR/ARGPARSE_INCLUDE_DIR."
        )
    endif()
endif()

# -----------------------------------------------------------------------------
# Sources and core target
# -----------------------------------------------------------------------------
file(GLOB_RECURSE FEMASTER_CPP_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
list(REMOVE_ITEM FEMASTER_CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

set(FEMASTER_SOURCES ${FEMASTER_CPP_SOURCES})
if(FEMASTER_ENABLE_CUDA)
    file(GLOB FEMASTER_CUDA_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/cuda/*.cu"
    )
    list(APPEND FEMASTER_SOURCES ${FEMASTER_CUDA_SOURCES})
endif()

add_library(femaster_core STATIC ${FEMASTER_SOURCES})
target_include_directories(femaster_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

if(Eigen3_FOUND)
    target_link_libraries(femaster_core PUBLIC Eigen3::Eigen)
else()
    target_include_directories(femaster_core PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()

if(Spectra_FOUND)
    target_link_libraries(femaster_core PUBLIC Spectra::Spectra)
else()
    target_include_directories(femaster_core PUBLIC "${SPECTRA_INCLUDE_DIR}")
endif()

if(argparse_FOUND)
    target_link_libraries(femaster_core PUBLIC argparse::argparse)
else()
    target_include_directories(femaster_core PUBLIC "${ARGPARSE_INCLUDE_DIR}")
endif()

target_compile_definitions(femaster_core
    PUBLIC
        $<$<BOOL:${FEMASTER_DOUBLE_PRECISION}>:DOUBLE_PRECISION>
        $<$<BOOL:${FEMASTER_DOUBLE_PRECISION}>:CUDA_DOUBLE_PRECISION>
        $<$<BOOL:${FEMASTER_FAST_EIGEN_COMPILE}>:EIGEN_UNROLLING_LIMIT=50>
        $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>:EIGEN_NO_DEBUG>
)

if(FEMASTER_TIME_REPORT)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
endif()

set(FEMASTER_MSVC_MP_JOBS "8" CACHE STRING "MSVC /MP parallel compiler processes")

if(MSVC)
    target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/MP${FEMASTER_MSVC_MP_JOBS}>
            $<$<COMPILE_LANGUAGE:CXX>:/bigobj>
            $<$<COMPILE_LANGUAGE:CXX>:/permissive->
    )
else()
    target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>
    )
endif()

# -----------------------------------------------------------------------------
# Intel oneMKL
# -----------------------------------------------------------------------------
set(FEMASTER_MKL_LINK "disabled")
set(FEMASTER_MKL_THREADING "disabled")
unset(FEMASTER_IOMP5_DLL CACHE)
unset(FEMASTER_IOMP5_DLL)

if(FEMASTER_ENABLE_MKL)
    set(MKL_INTERFACE "lp64")

    if(FEMASTER_MKL_STATIC)
        set(MKL_LINK "static")
    else()
        set(MKL_LINK "dynamic")
    endif()

    if(FEMASTER_MKL_SEQUENTIAL)
        set(MKL_THREADING "sequential")
    else()
        set(MKL_THREADING "intel_thread")
    endif()

    # Intel's package otherwise chooses libiomp5.so even for static MKL on
    # Linux. Select the static OpenMP runtime to match a distributable build.
    if(FEMASTER_MKL_STATIC AND NOT FEMASTER_MKL_SEQUENTIAL
       AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
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
                "Static threaded oneMKL requires libiomp5.a. Set "
                "FEMASTER_IOMP5_STATIC, enable FEMASTER_MKL_SEQUENTIAL, or "
                "disable FEMASTER_MKL_STATIC."
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
            message(FATAL_ERROR
                "Threaded oneMKL on Windows requires libiomp5md.dll, but it "
                "was not found in the oneAPI installation."
            )
        endif()
    endif()
endif()

# -----------------------------------------------------------------------------
# OpenMP
# -----------------------------------------------------------------------------
set(FEMASTER_OPENMP_RUNTIME "disabled")
if(FEMASTER_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_compile_definitions(femaster_core PUBLIC EIGEN_DONT_PARALLELIZE)

    if(FEMASTER_ENABLE_MKL AND NOT FEMASTER_MKL_SEQUENTIAL)
        if(MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            message(FATAL_ERROR
                "MSVC OpenMP and threaded oneMKL use different runtimes. "
                "Disable FEMASTER_ENABLE_OPENMP, select sequential MKL, or "
                "use the IntelLLVM compiler."
            )
        endif()

        # Compile FEMaster's OpenMP loops but let MKL::MKL provide Intel's
        # OpenMP runtime. This avoids linking libgomp and libiomp together.
        get_target_property(FEMASTER_OPENMP_COMPILE_OPTIONS
            OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS
        )
        if(FEMASTER_OPENMP_COMPILE_OPTIONS)
            target_compile_options(femaster_core PUBLIC
                ${FEMASTER_OPENMP_COMPILE_OPTIONS}
            )
        endif()
        set(FEMASTER_OPENMP_RUNTIME "Intel OpenMP via oneMKL")
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
    target_link_libraries(femaster_core
        PUBLIC
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

# -----------------------------------------------------------------------------
# Compiler options
# -----------------------------------------------------------------------------
if(MSVC)
    target_compile_options(femaster_core PRIVATE /bigobj /permissive-)
else()
    target_compile_options(femaster_core PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>
    )
endif()

if(FEMASTER_TIME_REPORT)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-ftime-report>
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(femaster_core PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-ftime-trace>
        )
    elseif(MSVC)
        target_compile_options(femaster_core PRIVATE /Bt+)
    endif()
endif()

if(MSVC AND FEMASTER_STATIC_RUNTIME)
    set_property(TARGET femaster_core PROPERTY MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>"
    )
endif()

# -----------------------------------------------------------------------------
# Executable and runtime packaging
# -----------------------------------------------------------------------------
add_executable(FEMaster src/main.cpp)
target_link_libraries(FEMaster PRIVATE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,femaster_core>"
)

if(MSVC AND FEMASTER_STATIC_RUNTIME)
    set_property(TARGET FEMaster PROPERTY MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>"
    )
elseif(FEMASTER_STATIC_RUNTIME AND NOT APPLE
       AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_options(FEMaster PRIVATE -static-libgcc -static-libstdc++)
endif()

if(FEMASTER_FULLY_STATIC)
    target_link_options(FEMaster PRIVATE -static)
endif()

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
            "${FEMASTER_IOMP5_DLL}"
            "$<TARGET_FILE_DIR:FEMaster>"
        COMMENT "Copying Intel OpenMP runtime"
    )
endif()

if(WIN32 AND FEMASTER_ENABLE_CUDSS)
    add_custom_command(TARGET FEMaster POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:cudss>"
            "$<TARGET_FILE_DIR:FEMaster>"
        COMMENT "Copying cuDSS runtime"
    )
endif()

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
if(FEMASTER_BUILD_TESTS)
    enable_testing()
    find_package(GTest CONFIG REQUIRED)

    file(GLOB TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp"
    )

    add_executable(FEMaster_test ${TEST_SOURCES})
    target_link_libraries(FEMaster_test
        PRIVATE
            "$<LINK_LIBRARY:WHOLE_ARCHIVE,femaster_core>"
            GTest::gtest
            GTest::gtest_main
    )
    if(MSVC AND FEMASTER_STATIC_RUNTIME)
        set_property(TARGET FEMaster_test PROPERTY MSVC_RUNTIME_LIBRARY
            "MultiThreaded$<$<CONFIG:Debug>:Debug>"
        )
    endif()

    include(GoogleTest)
    gtest_discover_tests(FEMaster_test)
endif()

# -----------------------------------------------------------------------------
# Human-readable configuration summary (also available as femaster-info target)
# -----------------------------------------------------------------------------
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
message(STATUS "  Double precision     : ${FEMASTER_DOUBLE_PRECISION}")
message(STATUS "  OpenMP loops         : ${FEMASTER_ENABLE_OPENMP} (${FEMASTER_OPENMP_RUNTIME})")
message(STATUS "  oneMKL               : ${FEMASTER_ENABLE_MKL}")
message(STATUS "    link/threading     : ${FEMASTER_MKL_LINK} / ${FEMASTER_MKL_THREADING}")
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
    COMMAND ${CMAKE_COMMAND} -E echo "OpenMP: ${FEMASTER_ENABLE_OPENMP} (${FEMASTER_OPENMP_RUNTIME})"
    COMMAND ${CMAKE_COMMAND} -E echo "oneMKL: ${FEMASTER_ENABLE_MKL} (${FEMASTER_MKL_LINK}, ${FEMASTER_MKL_THREADING})"
    COMMAND ${CMAKE_COMMAND} -E echo "CUDA: ${FEMASTER_ENABLE_CUDA} (${FEMASTER_CUDA_VERSION})"
    COMMAND ${CMAKE_COMMAND} -E echo "cuDSS: ${FEMASTER_ENABLE_CUDSS} (${FEMASTER_CUDSS_VERSION})"
    VERBATIM
)
