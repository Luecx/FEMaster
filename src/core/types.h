#pragma once

#include <cstdint>

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <memory>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

typedef uint64_t Time;
typedef  int32_t ID;
typedef uint8_t Dim;
typedef uint64_t Index;

// use double for all normal computations and cast to solver precision during solving
typedef double Precision;

#ifdef CUDA_DOUBLE_PRECISION
typedef double CudaPrecision;
#else
typedef float CudaPrecision;
#endif

// eigen redefinitions
template<size_t A, size_t B>
using StaticMatrix        = Eigen::Matrix<Precision, A             , B             >;
template<size_t B>
using SemiStaticMatrix    = Eigen::Matrix<Precision, Eigen::Dynamic, B             >;
using DynamicMatrix       = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
using MapMatrix           = Eigen::Map<DynamicMatrix>;
using SparseMatrix        = Eigen::SparseMatrix<Precision, Eigen::ColMajor>;
using Triplet             = Eigen::Triplet<Precision>;
using TripletList         = std::vector<Triplet>;

template<size_t A>
using StaticVector  = Eigen::Matrix<Precision, A, 1>;
using DynamicVector = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;

using Vec2          = StaticVector<2>;
using Vec3          = StaticVector<3>;
using Vec4          = StaticVector<4>;
using Vec6          = StaticVector<6>;

using IndexVector   = Eigen::Matrix<ID       , Eigen::Dynamic, 1>;
using IndexMatrix   = Eigen::Matrix<ID       , Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using BooleanVector = Eigen::Matrix<bool     , Eigen::Dynamic, 1>;
using BooleanMatrix = Eigen::Matrix<bool     , Eigen::Dynamic, Eigen::Dynamic>;

using NodeData      = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ElementData   = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Dofs          = Eigen::Matrix<bool     , 1             , 6             >;

using ElDofs       = Eigen::Matrix<bool, 1, 6>;     ///< Element DOFs: 1 row, 6 columns
using SystemDofs   = Eigen::Matrix<bool, Eigen::Dynamic, 6, Eigen::RowMajor>; ///< System DOFs: N rows, 6 columns
using SystemDofIds = Eigen::Matrix<int, Eigen::Dynamic, 6, Eigen::RowMajor>; ///< System DOF IDs: N rows, 6 columns

