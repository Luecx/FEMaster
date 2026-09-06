cmake_minimum_required(VERSION 3.24)

project(FRTShellTangentTest LANGUAGES CXX)

set(FEMASTER_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")

# Build only FEMaster's core library. The repository-wide test executable is
# intentionally disabled because unrelated legacy tests may not compile on the
# current branch.
set(FEMASTER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(FEMASTER_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
set(FEMASTER_ENABLE_CUDSS OFF CACHE BOOL "" FORCE)

add_subdirectory("${FEMASTER_ROOT}" "${CMAKE_BINARY_DIR}/femaster")

find_package(GTest CONFIG REQUIRED)

add_executable(FRTShellTangentTest
    "${FEMASTER_ROOT}/tests/test_shell_frt_tangent.cpp"
)

target_link_libraries(FRTShellTangentTest
    PRIVATE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,femaster_core>"
    GTest::gtest
    GTest::gtest_main
)

if(MSVC AND FEMASTER_STATIC_RUNTIME)
    set_property(TARGET FRTShellTangentTest PROPERTY MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>"
    )
endif()
