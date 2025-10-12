
cmake_minimum_required(VERSION 3.16)
project(FEMReaderV2 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(example examples/example_api.cpp)
target_include_directories(example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
