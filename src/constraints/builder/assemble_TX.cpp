/******************************************************************************
 * @file assemble_TX.cpp
 * @brief Assemble sparse T and X from X_cols and partition info.
 ******************************************************************************/

#include "assemble_TX.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <unordered_map>

namespace fem { namespace constraint {

AssembleOutput assemble_T_and_X(const AssembleInput& in, const XCols& Xc)
{
    AssembleOutput out{};
    const int nm_use = (int)in.masters_loc.size();
    const int nm     = (int)in.masters_glob.size();

    // ------------------------- Build T (n x nm) ------------------------------
    {
        std::vector<Eigen::Triplet<Precision>> trips;
        size_t nnz_X_est = 0;
        for (int j = 0; j < nm_use; ++j) nnz_X_est += Xc.cols[j].size();
        trips.reserve((size_t)nm + nnz_X_est);

        // Map constrained master global id -> local QR-master index jloc
        std::unordered_map<int,int> masterGlobalToLocal;
        masterGlobalToLocal.reserve((size_t)nm_use * 2);
        for (int j = 0; j < nm_use; ++j) {
            const int cg = in.used[in.masters_loc[j]];
            masterGlobalToLocal[cg] = j;
        }

        int col_idx = 0;
        for (; col_idx < nm; ++col_idx) {
            const int col_global = in.masters_glob[col_idx];
            auto itFind = masterGlobalToLocal.find(col_global);
            if (itFind == masterGlobalToLocal.end()) break; // trailing are never-used masters

            const int jloc = itFind->second;

            // 1 on master row
            trips.emplace_back(col_global, col_idx, Precision(1));

            // add slave rows per X column
            const auto& xj = Xc.cols[jloc];
            for (const auto& ij : xj) {
                const int i_slave_local = ij.first;
                const int slave_global  = in.used[in.slaves_loc[i_slave_local]];
                trips.emplace_back(slave_global, col_idx, ij.second);
            }
        }

        // identities for never-used masters
        for (; col_idx < nm; ++col_idx) {
            const int g = in.masters_glob[col_idx];
            trips.emplace_back(g, col_idx, Precision(1));
        }

        out.T.resize(in.n, nm);
        out.T.setFromTriplets(trips.begin(), trips.end());
        out.T.makeCompressed();
    }

    // ------------------------- Build X (r x nm) ------------------------------
    {
        out.X.resize(in.r, nm);
        if (in.r > 0 && nm_use > 0) {
            std::vector<Eigen::Triplet<Precision>> xTrips;
            size_t nnz_X = 0;
            for (int j = 0; j < nm_use; ++j) nnz_X += Xc.cols[j].size();
            xTrips.reserve(nnz_X);

            for (int j = 0; j < nm_use; ++j) {
                const auto& xj = Xc.cols[j];
                for (const auto& ij : xj) {
                    xTrips.emplace_back(ij.first, j, ij.second);
                }
            }
            out.X.setFromTriplets(xTrips.begin(), xTrips.end());
            out.X.makeCompressed();
        }
    }

    return out;
}

}} // namespace fem::constraint
