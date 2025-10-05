/******************************************************************************
 * @file partition.cpp
 * @brief Splits permutation data into slave and master DOF sets.
 *
 * The partitioning step maps local QR column indices back to global DOFs while
 * filtering out fixed entries.
 *
 * @see src/constraints/builder/partition.h
 * @see src/constraints/builder/preprocess.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "partition.h"

#include "../../core/logging.h"

#include <algorithm>

namespace fem {
namespace constraint {

/******************************************************************************
 * @copydoc partition_and_map
 ******************************************************************************/
PartitionOutput partition_and_map(const PartitionInput& input, const std::vector<char>& is_used) {
    PartitionOutput output;

    const auto& indices = input.P.indices();

    output.slaves_loc.reserve(input.r);
    for (int i = 0; i < input.r; ++i) {
        output.slaves_loc.push_back(indices(i));
    }

    const int n_use = static_cast<int>(input.used->size());
    output.nm_use = n_use - input.r;
    output.masters_loc.reserve(std::max(0, output.nm_use));
    for (int j = 0; j < output.nm_use; ++j) {
        output.masters_loc.push_back(indices(input.r + j));
    }

    output.slaves_glob.reserve(input.r);
    for (int i = 0; i < input.r; ++i) {
        const int cg = input.used->at(output.slaves_loc[i]);
        if (!input.is_fixed_col->at(cg)) {
            output.slaves_glob.push_back(cg);
        }
    }

    output.masters_glob.reserve(input.n);
    for (int j = 0; j < output.nm_use; ++j) {
        const int cg = input.used->at(output.masters_loc[j]);
        if (!input.is_fixed_col->at(cg)) {
            output.masters_glob.push_back(cg);
        }
    }

    for (int g = 0; g < input.n; ++g) {
        if (!is_used[g] && !input.is_fixed_col->at(g)) {
            output.masters_glob.push_back(g);
        }
    }

    return output;
}

} // namespace constraint
} // namespace fem
