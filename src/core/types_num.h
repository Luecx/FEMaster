//
// Created by f_eggers on 20.12.2024.
//

#ifndef TYPES_NUM_H
#define TYPES_NUM_H

#include <cstdint>


namespace fem {

// Type aliases for basic index and ID types
typedef uint64_t Time;  ///< Time type used for representing time in microseconds
typedef int32_t ID;     ///< ID type for identifying nodes, elements, etc.
typedef uint8_t Dim;    ///< Dimension type (e.g., 2D or 3D)
typedef uint64_t Index; ///< Index type for large indexing in arrays

/******************************************************************************
* @brief Precision typedefs for floating-point operations.
*
* The precision type for normal CPU computations is set by the DOUBLE_PRECISION flag.
* If DOUBLE_PRECISION is defined, Precision is double; otherwise, it is float. For
* CUDA-based computations, CudaPrecision follows a similar pattern.
******************************************************************************/
#ifdef DOUBLE_PRECISION
typedef double Precision;  ///< Precision for CPU computations
#else
typedef float Precision;   ///< Single precision for CPU computations
#endif

#ifdef CUDA_DOUBLE_PRECISION
typedef double CudaPrecision;  ///< Double precision for CUDA computations
#else
typedef float CudaPrecision;   ///< Single precision for CUDA computations
#endif

}  // namespace fem

#endif //TYPES_NUM_H
