#include "cb_partition.h"

namespace fem { namespace constraint {

PartitionOut partition_by_P(const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int>& P,
                            int r, int nm_use)
{
    PartitionOut out;
    out.slaves_loc.reserve(r);
    out.masters_loc.reserve(nm_use);

    const auto& p = P.indices();  // P * e_i = e_{p(i)}

    for (int i = 0; i < r; ++i)
        out.slaves_loc.push_back(p(i));

    for (int j = 0; j < nm_use; ++j)
        out.masters_loc.push_back(p(r + j));

    return out;
}

}} // namespace fem::constraint
