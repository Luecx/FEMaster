#pragma once
/**
 * @file cb_partition.h
 * @brief Extract slave/master local column partitions from the QR column permutation.
 *
 * After SparseQR on C_use, we use the permutation P such that C_use * P = Q * R.
 * If R has numerical rank r, then the first r permuted columns correspond to
 * “slave” columns and the next nm_use columns correspond to “master” columns.
 *
 * This helper gathers those sets as *local indices* into the `used` array
 * (i.e., indices w.r.t. C_use’s columns).
 */

#include <Eigen/Sparse>
#include <vector>

namespace fem { namespace constraint {

struct PartitionOut {
    std::vector<int> slaves_loc;   ///< size r; local indices in [0..n_use-1]
    std::vector<int> masters_loc;  ///< size nm_use; local indices in [0..n_use-1]
};

/**
 * @brief Read the column permutation to produce local slave/master partitions.
 *
 * @param P       Column permutation from SparseQR (so that C_use * P = Q * R).
 * @param r       Numerical rank (number of slave columns).
 * @param nm_use  Number of master columns (n_use - r).
 * @return PartitionOut with local indices for slaves and masters.
 */
PartitionOut partition_by_P(const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int>& P,
                            int r, int nm_use);

}} // namespace fem::constraint
