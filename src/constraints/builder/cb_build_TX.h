#pragma once
/**
 * @file cb_build_TX.h
 * @brief Assemble the reduced null-space transformer T and the slave-coupling X.
 *
 * Given the QR partition of the constrained columns, this builds:
 *   - T_full (n × nm_total): maps reduced coordinates q to full u via u = u_p + T q
 *   - X      (r × nm_total): top-left block equals Xd = -R11^{-1} R12 (slave ← master)
 *
 * Inputs:
 *   n              : total number of DOFs (global columns of C)
 *   r              : numerical rank of C_use (number of slave DOFs from QR)
 *   used           : global column indices retained in C_use after compression
 *   slaves_loc     : length r; indices into `used` that are slave columns (local QR order)
 *   masters_loc    : length nm_use; indices into `used` that are master columns (local QR order)
 *   is_used        : size n; marks columns that appear in C_use
 *   is_fixed_col   : size n; columns fixed by single-NNZ rows (become part of u_p, not masters)
 *   X_cols         : sparse columns of Xd (= -R11^{-1}R12), one per QR master (local)
 *
 * Outputs:
 *   TXOut {
 *     T_full       : n × nm_total (nm_total = #QR masters + #never-used masters), sparse
 *     X            : r × nm_total (left block is Xd, trailing never-used masters are zero)
 *     masters_glob : global indices of master DOFs (order matches columns of T_full/X)
 *     slaves_glob  : global indices of slave DOFs (order matches rows 0..r-1 of X)
 *   }
 *
 * Notes:
 *   - “never-used masters” are the global columns that never appear in any kept row of C
 *     (i.e., !is_used[g]) and are also not fixed; they map as identity columns in T.
 *   - Fixed columns never appear as masters or slaves; they contribute only to u_p.
 */

#include <Eigen/Sparse>
#include <unordered_map>
#include <vector>
#include "../../core/types_eig.h"
#include "cb_build_x_sparse.h"    // for XCols

namespace fem { namespace constraint {

struct TXOut {
    /// n × nm_total. Each column j has a 1 at its master row plus X entries on slave rows.
    SparseMatrix T_full;

    /// r × nm_total. Left block equals Xd (one column per QR master). Zero for never-used masters.
    SparseMatrix X;

    /// Global indices of master DOFs (columns of T_full/X).
    std::vector<Index> masters_glob;

    /// Global indices of slave DOFs (rows 0..r-1 of X refer to these).
    std::vector<Index> slaves_glob;
};

/**
 * @brief Build T_full and X from QR partition and sparse Xd columns.
 *
 * @param n             Total number of DOFs (global dimension).
 * @param r             Rank (number of slave DOFs).
 * @param used          Global indices of columns retained in C_use (size = n_use).
 * @param slaves_loc    Local QR indices (into `used`) of slave columns, size r.
 * @param masters_loc   Local QR indices (into `used`) of master columns, size nm_use.
 * @param is_used       Flags of size n marking columns present in C_use.
 * @param is_fixed_col  Flags of size n marking columns fixed by simple constraints.
 * @param X_cols        Sparse columns for Xd, length nm_use; X_cols[j] holds (i,val) pairs with i in [0..r-1].
 * @return TXOut        Assembled T_full, X, and the global master/slave index lists.
 */
TXOut build_T_and_X(int n,
                    int r,
                    const std::vector<int>& used,
                    const std::vector<int>& slaves_loc,
                    const std::vector<int>& masters_loc,
                    const std::vector<char>& is_used,
                    const std::vector<char>& is_fixed_col,
                    const XCols& X_cols);

}} // namespace fem::constraint
