/**
 * @file types_eig.h
 * @brief Declares Eigen-based matrix and vector aliases used across the solver.
 *
 * Collects commonly used dense, sparse, mapped, and fixed-size Eigen types
 * using the configured numeric precision.
 *
 * @see src/core/types_num.h
 * @see src/core/types_cls.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#include "types_num.h"

namespace fem {

// ============================================================================
// Generic dense types
// ============================================================================

template<int Rows, int Cols>
using StaticMatrix     = Eigen::Matrix<Precision, Rows, Cols>;

template<int Cols>
using SemiStaticMatrix = Eigen::Matrix<Precision, Eigen::Dynamic, Cols>;

template<int Rows>
using StaticVector     = Eigen::Matrix<Precision, Rows, 1>;

using DynamicMatrix    = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
using DynamicVector    = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;
using RowMatrix        = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MapMatrix        = Eigen::Map<DynamicMatrix>;

// ============================================================================
// Sparse types
// ============================================================================

using SparseMatrix = Eigen::SparseMatrix<Precision, Eigen::ColMajor>;
using Triplet      = Eigen::Triplet<Precision>;
using TripletList  = std::vector<Triplet>;

// ============================================================================
// Fixed-size vectors
// ============================================================================

using Vec2 = StaticVector<2>;
using Vec3 = StaticVector<3>;
using Vec4 = StaticVector<4>;
using Vec5 = StaticVector<5>;
using Vec6 = StaticVector<6>;
using Vec8 = StaticVector<8>;

// ============================================================================
// Fixed-size matrices
// ============================================================================

using Mat2 = StaticMatrix<2, 2>;
using Mat3 = StaticMatrix<3, 3>;
using Mat4 = StaticMatrix<4, 4>;
using Mat5 = StaticMatrix<5, 5>;
using Mat6 = StaticMatrix<6, 6>;
using Mat8 = StaticMatrix<8, 8>;

using Mat23 = StaticMatrix<2, 3>;
using Mat24 = StaticMatrix<2, 4>;
using Mat32 = StaticMatrix<3, 2>;
using Mat34 = StaticMatrix<3, 4>;
using Mat42 = StaticMatrix<4, 2>;
using Mat43 = StaticMatrix<4, 3>;

// ============================================================================
// Index and boolean types
// ============================================================================

using IndexVector   = Eigen::Matrix<ID, Eigen::Dynamic, 1>;
using IndexMatrix   = Eigen::Matrix<ID, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using BooleanVector = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using BooleanMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

// ============================================================================
// Degree-of-freedom types
// ============================================================================

using Dofs         = Eigen::Matrix<bool, 1, 6>;
using ElDofs       = Dofs;
using SystemDofs   = Eigen::Matrix<bool, Eigen::Dynamic, 6, Eigen::RowMajor>;
using SystemDofIds = Eigen::Matrix<int, Eigen::Dynamic, 6, Eigen::RowMajor>;

// ============================================================================
// Integration field types
// ============================================================================

using ScalarField = std::function<Precision(const Vec3&)>;
using VecField    = std::function<Vec3(const Vec3&)>;
using TenField    = std::function<Mat3(const Vec3&)>;

} // namespace fem