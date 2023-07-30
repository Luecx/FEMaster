#pragma once

#include <cstdint>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigen>

typedef uint64_t Time;
typedef uint32_t ID;
typedef uint8_t Dim;

#ifdef DOUBLE_PRECISION
typedef double Precision;
#else
typedef float Precision;
#endif

// eigen redefinitions
template<size_t A, size_t B>
using StaticMatrix        = Eigen::Matrix<Precision, A             , B             >;
using DynamicMatrix       = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMatrix        = Eigen::SparseMatrix<Precision>;
using SparseMatrixBuilder = std::vector<Eigen::Triplet<Precision>>;
template<size_t A>
using StaticVector  = Eigen::Matrix<Precision, A, 1>;
using DynamicVector = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;

using IndexVector   = Eigen::Matrix<ID       , Eigen::Dynamic, 1>;
using IndexMatrix   = Eigen::Matrix<ID       , Eigen::Dynamic, Eigen::Dynamic>;

using BooleanVector = Eigen::Matrix<bool     , Eigen::Dynamic, 1>;
using BooleanMatrix = Eigen::Matrix<bool     , Eigen::Dynamic, Eigen::Dynamic>;

using NodeData      = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
using Dofs          = Eigen::Matrix<bool     , 6             , 1             >;

