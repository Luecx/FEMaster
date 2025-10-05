/**
 * @file assemble_TX.cpp
 * @brief Builds the `T` and `X` matrices from reduced column data.
 *
 * This step combines output from the QR factorisation with partitioning to
 * construct the null-space basis and explicit constraint matrix.
 *
 * @see src/constraints/builder/assemble_TX.h
 * @see src/constraints/builder/build_x.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "assemble_TX.h"

#include <unordered_map>

namespace fem {
namespace constraint {

/**
 * @copydoc assemble_T_and_X
 */
AssembleOutput assemble_T_and_X(const AssembleInput& in, const XCols& cols) {
    AssembleOutput out;
    const int nm_use = static_cast<int>(in.masters_loc.size());
    const int nm = static_cast<int>(in.masters_glob.size());

    std::vector<Eigen::Triplet<Precision>> trips;
    size_t nnz_estimate = 0;
    for (int j = 0; j < nm_use; ++j) {
        nnz_estimate += cols.cols[j].size();
    }
    trips.reserve(static_cast<std::size_t>(nm) + nnz_estimate);

    std::unordered_map<int, int> master_global_to_local;
    master_global_to_local.reserve(static_cast<std::size_t>(nm_use) * 2);
    for (int j = 0; j < nm_use; ++j) {
        const int cg = in.used[in.masters_loc[j]];
        master_global_to_local[cg] = j;
    }

    int column_index = 0;
    for (; column_index < nm; ++column_index) {
        const int column_global = in.masters_glob[column_index];
        auto found = master_global_to_local.find(column_global);
        if (found == master_global_to_local.end()) {
            break;
        }

        const int local_index = found->second;
        trips.emplace_back(column_global, column_index, Precision(1));

        const auto& column_entries = cols.cols[local_index];
        for (const auto& entry : column_entries) {
            const int slave_local = entry.first;
            const int slave_global = in.used[in.slaves_loc[slave_local]];
            trips.emplace_back(slave_global, column_index, entry.second);
        }
    }

    for (; column_index < nm; ++column_index) {
        const int g = in.masters_glob[column_index];
        trips.emplace_back(g, column_index, Precision(1));
    }

    out.T.resize(in.n, nm);
    out.T.setFromTriplets(trips.begin(), trips.end());
    out.T.makeCompressed();

    out.X.resize(in.r, nm);
    if (in.r > 0 && nm_use > 0) {
        std::vector<Eigen::Triplet<Precision>> x_trips;
        size_t nnz_X = 0;
        for (int j = 0; j < nm_use; ++j) {
            nnz_X += cols.cols[j].size();
        }
        x_trips.reserve(nnz_X);

        for (int j = 0; j < nm_use; ++j) {
            const auto& column_entries = cols.cols[j];
            for (const auto& entry : column_entries) {
                x_trips.emplace_back(entry.first, j, entry.second);
            }
        }
        out.X.setFromTriplets(x_trips.begin(), x_trips.end());
        out.X.makeCompressed();
    }

    return out;
}

} // namespace constraint
} // namespace fem
