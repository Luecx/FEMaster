/******************************************************************************
* @file partition.h
 * @brief Split permutation P into local slaves/masters and map to globals.
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"
#include <vector>

namespace fem { namespace constraint {

struct PartitionInput {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
    int r;                         ///< rank
    const std::vector<int>* used;  ///< local→global map for C_use columns
    const std::vector<char>* is_fixed_col; ///< global mask
    int n;                         ///< total columns
};

struct PartitionOutput {
    std::vector<int> slaves_loc;     ///< size r (local indices in used)
    std::vector<int> masters_loc;    ///< size nm_use
    std::vector<Index> slaves_glob;  ///< filtered (exclude fixed)
    std::vector<Index> masters_glob; ///< constrained masters + never-used
    int nm_use = 0;
};

PartitionOutput partition_and_map(const PartitionInput& in, const std::vector<char>& is_used);

}} // namespace
