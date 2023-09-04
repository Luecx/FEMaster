cmake_minimum_required(VERSION 3.17)
project(FEMaster VERSION 1.0 LANGUAGES CUDA CXX)

find_package(CUDA REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CUDA_STANDARD 17)

file(GLOB_RECURSE CPP_SRCS src/*.cpp)
file(GLOB_RECURSE CU_SRCS src/*.cu)

## CPU target
#add_executable(FEMaster_CPU ${CPP_SRCS})
#set_target_properties(FEMaster_CPU PROPERTIES
#        CUDA_SEPARABLE_COMPILATION ON
#        )
#
#target_compile_options(FEMaster_CPU PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)
#
#target_link_libraries(FEMaster_CPU PRIVATE
#        ${CMAKE_THREAD_LIBS_INIT}
#        )

# GPU target
add_executable(FEMaster_GPU ${CPP_SRCS} ${CU_SRCS})

set_property(TARGET FEMaster_GPU PROPERTY CUDA_ARCHITECTURES OFF)
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")

set_target_properties(FEMaster_GPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -openmp ")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 ")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -use_fast_math")

target_compile_definitions(FEMaster_GPU PRIVATE SUPPORT_GPU)

target_compile_options(FEMaster_GPU PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)

include_directories("include")

target_link_libraries(FEMaster_GPU ${CUDA_LIBRARIES})
target_link_libraries(FEMaster_GPU ${CUDA_CUBLAS_LIBRARIES})
target_link_libraries(FEMaster_GPU ${CUDA_cusparse_LIBRARY})
target_link_libraries(FEMaster_GPU ${CMAKE_THREAD_LIBS_INIT})

#target_link_libraries(FEMaster_GPU
#        ${CUDA_LIBRARIES}
#        ${CUDA_CUBLAS_LIBRARIES}
#        ${CUDA_cusparse_LIBRARY}
#        ${CMAKE_THREAD_LIBS_INIT}
#        ${CUDA_cusolver_LIBRARY}
#        cusolver
#        )
