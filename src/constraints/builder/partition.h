/******************************************************************************
 * @file partition.h
 * @brief Declares utilities to split permutation data into slave/master sets.
 *
 * Partitioning separates constrained (slave) and unconstrained (master) DOFs
 * and maps them back to global indices.
 *
 * @see src/constraints/builder/partition.cpp
 * @see src/constraints/builder/preprocess.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../../core/types_eig.h"

#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @struct PartitionInput
 * @brief Data required to partition QR column permutations.
 ******************************************************************************/
struct PartitionInput {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P; ///< Column permutation from QR.
    int r = 0;                                ///< Numerical rank.
    const std::vector<int>* used = nullptr;   ///< Mapping from local to global columns.
    const std::vector<char>* is_fixed_col = nullptr; ///< Global mask of fixed DOFs.
    int n = 0;                                ///< Total number of columns.
};

/******************************************************************************
 * @struct PartitionOutput
 * @brief Partition results including local and global mappings.
 ******************************************************************************/
struct PartitionOutput {
    std::vector<int> slaves_loc;     ///< Local indices of slave columns.
    std::vector<int> masters_loc;    ///< Local indices of master columns.
    std::vector<Index> slaves_glob;  ///< Global indices of slave columns.
    std::vector<Index> masters_glob; ///< Global indices of master columns (including unused).
    int nm_use = 0;                  ///< Number of masters participating in the QR factor.
};

/******************************************************************************
 * @brief Partitions the permutation into slave/master sets and maps to globals.
 *
 * @param input Partitioning input data.
 * @param is_used Flags indicating which global columns are part of the reduced system.
 * @return PartitionOutput Partitioned indices and mappings.
 ******************************************************************************/
PartitionOutput partition_and_map(const PartitionInput& input, const std::vector<char>& is_used);

} // namespace constraint
} // namespace fem
