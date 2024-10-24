/******************************************************************************
* @file types.h
* @brief Provides type definitions for matrices, vectors, and precision handling
* used in the FEMaster solver. This includes typedefs for Eigen types, precision
* settings, and indexing types, ensuring flexibility for double or single precision.
*
* The file defines matrix and vector types based on Eigen, as well as precision
* settings for CPU and CUDA-based computations. These types are used throughout
* the FEMaster solver to maintain consistency and readability.
*
* @note The precision type can be set to either float or double, depending on the
* definition of the DOUBLE_PRECISION flag during compilation.
*
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#pragma once  // Ensures this file is only included once during compilation

#include <cstdint>

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL  // Use Intel MKL for Eigen if available
#endif

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <memory>  // For std::shared_ptr

#ifdef USE_MKL
#include <Eigen/PardisoSupport>  // Use Pardiso solver when MKL is enabled
#endif

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

/******************************************************************************
* @brief Type aliases for Eigen matrix and vector definitions.
*
* These type aliases wrap Eigen matrix and vector types using the defined Precision
* type. StaticMatrix, SemiStaticMatrix, and DynamicMatrix are used for matrices with
* varying dimensions, while MapMatrix is used to map existing memory into an Eigen matrix.
******************************************************************************/
template<size_t A, size_t B>
using StaticMatrix = Eigen::Matrix<Precision, A, B>;                 ///< Fixed-size matrix
template<size_t B>
using SemiStaticMatrix = Eigen::Matrix<Precision, Eigen::Dynamic, B>; ///< Semi-static matrix with dynamic rows
using DynamicMatrix = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>; ///< Fully dynamic matrix
using MapMatrix = Eigen::Map<DynamicMatrix>; ///< Eigen matrix mapped to external memory
using SparseMatrix = Eigen::SparseMatrix<Precision, Eigen::ColMajor>; ///< Sparse matrix for FEM assembly
using Triplet = Eigen::Triplet<Precision>; ///< Triplet structure for sparse matrix construction
using TripletList = std::vector<Triplet>;  ///< List of triplets for sparse matrix assembly

/******************************************************************************
* @brief Type aliases for Eigen vectors of varying sizes.
*
* StaticVector defines fixed-size vectors, while DynamicVector is used for vectors
* with dynamic sizes. Vec2, Vec3, Vec4, and Vec6 are common fixed-size vectors.
******************************************************************************/
template<size_t A>
using StaticVector = Eigen::Matrix<Precision, A, 1>; ///< Fixed-size vector
using DynamicVector = Eigen::Matrix<Precision, Eigen::Dynamic, 1>; ///< Dynamic-size vector

using Vec2 = StaticVector<2>; ///< 2D vector
using Vec3 = StaticVector<3>; ///< 3D vector
using Vec4 = StaticVector<4>; ///< 4D vector
using Vec6 = StaticVector<6>; ///< 6D vector

/******************************************************************************
* @brief Index and Boolean matrix types for FEM system DOFs and IDs.
*
* These types are used to store node indices, system DOF IDs, and Boolean matrices
* for active DOFs in FEM systems. SystemDofs and SystemDofIds are defined as
* dynamic-size matrices with 6 columns, representing the degrees of freedom per node.
******************************************************************************/
using IndexVector = Eigen::Matrix<ID, Eigen::Dynamic, 1>; ///< Vector of element/node IDs
using IndexMatrix = Eigen::Matrix<ID, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; ///< Index matrix in row-major order

using BooleanVector = Eigen::Matrix<bool, Eigen::Dynamic, 1>; ///< Dynamic-size Boolean vector
using BooleanMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>; ///< Dynamic-size Boolean matrix

/******************************************************************************
* @brief Data types for storing node and element data.
*
* NodeData and ElementData are row-major matrices used to store data related to
* nodes and elements in the FEM system.
******************************************************************************/
using NodeData = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; ///< Node-related data matrix
using ElementData = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; ///< Element-related data matrix

/******************************************************************************
* @brief DOF-related matrix definitions for FEM systems.
*
* These matrices store information about the degrees of freedom (DOFs) for elements
* and systems in the FEM model. ElDofs represents the DOFs for a single element,
* while SystemDofs and SystemDofIds handle system-wide DOFs.
******************************************************************************/
using Dofs = Eigen::Matrix<bool, 1, 6>; ///< Element DOFs (1 row, 6 columns)
using ElDofs = Eigen::Matrix<bool, 1, 6>; ///< Element DOFs for a single element
using SystemDofs = Eigen::Matrix<bool, Eigen::Dynamic, 6, Eigen::RowMajor>; ///< System DOFs for all elements/nodes
using SystemDofIds = Eigen::Matrix<int, Eigen::Dynamic, 6, Eigen::RowMajor>; ///< System DOF IDs for mapping local to global DOFs

} // namespace fem