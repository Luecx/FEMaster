/******************************************************************************
* @file partition.cpp
 * @brief Split permutation P into local slaves/masters and map to global indices.
 ******************************************************************************/

#include "partition.h"

#include "../../core/logging.h"
#include <algorithm>

namespace fem { namespace constraint {

PartitionOutput partition_and_map(const PartitionInput& in, const std::vector<char>& is_used)
{
    PartitionOutput out{};

    const auto& p = in.P.indices();

    // Local indices according to P: first r are slaves, remaining are masters (within C_use)
    out.slaves_loc.reserve(in.r);
    for (int i = 0; i < in.r; ++i) out.slaves_loc.push_back(p(i));

    const int n_use = (int)in.used->size();
    out.nm_use = n_use - in.r;
    out.masters_loc.reserve(std::max(0, out.nm_use));
    for (int j = 0; j < out.nm_use; ++j) out.masters_loc.push_back(p(in.r + j));

    // Map to global and filter fixed columns
    out.slaves_glob.reserve(in.r);
    for (int i = 0; i < in.r; ++i) {
        const int cg = in.used->at(out.slaves_loc[i]);
        if (!in.is_fixed_col->at(cg)) out.slaves_glob.push_back(cg);
    }

    // masters_glob: constrained masters + never-used (identity), both excluding fixed
    out.masters_glob.reserve(in.n);

    // constrained masters (from QR)
    for (int j = 0; j < out.nm_use; ++j) {
        const int cg = in.used->at(out.masters_loc[j]);
        if (!in.is_fixed_col->at(cg)) out.masters_glob.push_back(cg);
    }

    // append never-used columns (not in used, not fixed)
    for (int g = 0; g < in.n; ++g) {
        if (!is_used[g] && !in.is_fixed_col->at(g)) out.masters_glob.push_back(g);
    }

    return out;
}

}} // namespace fem::constraint
