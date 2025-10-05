/******************************************************************************
 * @file types_eig.h
 * @brief Declares Eigen-based matrix and vector aliases used across the solver.
 *
 * Collects commonly used dense, sparse, and map types tied to the configured
 * numeric precision.
 *
 * @see src/core/types_num.h
 * @see src/core/types_cls.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include <cstdint>
#include <memory>
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

template<std::size_t Rows, std::size_t Cols>
using StaticMatrix = Eigen::Matrix<Precision, Rows, Cols>;

template<std::size_t Cols>
using SemiStaticMatrix = Eigen::Matrix<Precision, Eigen::Dynamic, Cols>;

using DynamicMatrix = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
using MapMatrix = Eigen::Map<DynamicMatrix>;
using SparseMatrix = Eigen::SparseMatrix<Precision, Eigen::ColMajor>;
using Triplet = Eigen::Triplet<Precision>;
using TripletList = std::vector<Triplet>;

template<std::size_t Rows>
using StaticVector = Eigen::Matrix<Precision, Rows, 1>;

using DynamicVector = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;
using Vec2 = StaticVector<2>;
using Vec3 = StaticVector<3>;
using Vec4 = StaticVector<4>;
using Vec6 = StaticVector<6>;

using Mat2 = StaticMatrix<2, 2>;
using Mat3 = StaticMatrix<3, 3>;
using Mat4 = StaticMatrix<4, 4>;
using Mat6 = StaticMatrix<6, 6>;
using Mat42 = StaticMatrix<4, 2>;
using Mat24 = StaticMatrix<2, 4>;
using Mat32 = StaticMatrix<3, 2>;
using Mat23 = StaticMatrix<2, 3>;
using Mat43 = StaticMatrix<4, 3>;
using Mat34 = StaticMatrix<3, 4>;

using IndexVector = Eigen::Matrix<ID, Eigen::Dynamic, 1>;
using IndexMatrix = Eigen::Matrix<ID, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using BooleanVector = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using BooleanMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

using NodeData = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ElementData = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using IPData = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<std::size_t Cols>
using NodeDataM = Eigen::Matrix<Precision, Eigen::Dynamic, Cols, Eigen::RowMajor>;

template<std::size_t Cols>
using ElementDataM = Eigen::Matrix<Precision, Eigen::Dynamic, Cols, Eigen::RowMajor>;

template<std::size_t Cols>
using IPDataM = Eigen::Matrix<Precision, Eigen::Dynamic, Cols, Eigen::RowMajor>;

using Dofs = Eigen::Matrix<bool, 1, 6>;
using ElDofs = Eigen::Matrix<bool, 1, 6>;
using SystemDofs = Eigen::Matrix<bool, Eigen::Dynamic, 6, Eigen::RowMajor>;
using SystemDofIds = Eigen::Matrix<int, Eigen::Dynamic, 6, Eigen::RowMajor>;

} // namespace fem
