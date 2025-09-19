#include "cb_build_TX.h"

namespace fem { namespace constraint {

TXOut build_T_and_X(int n,
                    int r,
                    const std::vector<int>& used,
                    const std::vector<int>& slaves_loc,
                    const std::vector<int>& masters_loc,
                    const std::vector<char>& is_used,
                    const std::vector<char>& is_fixed_col,
                    const XCols& X_cols)
{
    const int nm_use = static_cast<int>(masters_loc.size());

    // -------------------------------------------------------------------------
    // 1) Build master and slave global index lists
    //    - masters_glob: QR masters first (excluding fixed), then “never-used” masters (excluding fixed)
    //    - slaves_glob : QR slaves (excluding fixed, which shouldn't appear anyway)
    // -------------------------------------------------------------------------
    std::vector<Index> masters_glob;
    masters_glob.reserve(n);

    // QR masters (map local QR indices -> global col in `used`)
    for (int j = 0; j < nm_use; ++j) {
        const int col_local  = masters_loc[j];
        const int col_global = used[col_local];
        if (!is_fixed_col[col_global]) {
            masters_glob.push_back(col_global);
        }
    }

    // “Never-used” masters: columns not present in C_use and not fixed → identity columns in T
    for (int g = 0; g < n; ++g) {
        if (!is_used[g] && !is_fixed_col[g]) {
            masters_glob.push_back(g);
        }
    }

    std::vector<Index> slaves_glob;
    slaves_glob.reserve(r);
    for (int i = 0; i < r; ++i) {
        const int col_global = used[slaves_loc[i]];
        if (!is_fixed_col[col_global]) {
            slaves_glob.push_back(col_global);
        }
    }

    const int nm_total = static_cast<int>(masters_glob.size());

    // -------------------------------------------------------------------------
    // 2) Build T_full (n × nm_total), as triplets:
    //      - For each QR master column j: a 1 at its master row plus Xd entries on slave rows.
    //      - For each never-used master: a single identity entry at (g,g_col).
    //
    //    We need a quick map: global master → local QR master j (0..nm_use-1).
    //    For never-used masters this map misses, and we emit identity later.
    // -------------------------------------------------------------------------
    SparseMatrix T_full(n, nm_total);
    {
        std::vector<Eigen::Triplet<Precision>> t_trips;
        // Rough reserve: 1 per master + all Xd nnz
        std::size_t nnz_x_est = 0;
        for (int j = 0; j < nm_use; ++j) nnz_x_est += X_cols[j].size();
        t_trips.reserve(static_cast<std::size_t>(nm_total) + nnz_x_est);

        std::unordered_map<int, int> master_global_to_local;
        master_global_to_local.reserve(static_cast<std::size_t>(nm_use) * 2);
        for (int j = 0; j < nm_use; ++j) {
            master_global_to_local.emplace(used[masters_loc[j]], j);
        }

        // Emit QR master columns first (same order as we placed them into masters_glob)
        int col_out = 0;
        for (; col_out < nm_total; ++col_out) {
            const int g_master = static_cast<int>(masters_glob[col_out]);

            auto it = master_global_to_local.find(g_master);
            if (it == master_global_to_local.end()) break;  // “never-used” masters start here

            const int jloc = it->second;  // local QR-master index [0..nm_use-1]

            // 1 at the master row
            t_trips.emplace_back(g_master, col_out, Precision(1));

            // Add Xd entries on the corresponding slave rows
            const auto& xj = X_cols[jloc];  // vector of (i_slave_local, value)
            for (const auto& ij : xj) {
                const int i_slave_local = ij.first;  // 0..r-1
                const int g_slave       = used[slaves_loc[i_slave_local]];
                const Precision v       = ij.second;
                t_trips.emplace_back(g_slave, col_out, v);
            }
        }

        // Emit identities for trailing “never-used” masters
        for (; col_out < nm_total; ++col_out) {
            const int g = static_cast<int>(masters_glob[col_out]);
            t_trips.emplace_back(g, col_out, Precision(1));
        }

        T_full.setFromTriplets(t_trips.begin(), t_trips.end());
        T_full.makeCompressed();
    }

    // -------------------------------------------------------------------------
    // 3) Build X (r × nm_total):
    //      - Left block (first nm_use columns) = Xd from X_cols
    //      - Trailing columns (never-used masters) = zeros
    // -------------------------------------------------------------------------
    SparseMatrix X(r, nm_total);
    {
        std::vector<Eigen::Triplet<Precision>> x_trips;
        std::size_t nnz_x = 0;
        for (int j = 0; j < nm_use; ++j) nnz_x += X_cols[j].size();
        x_trips.reserve(nnz_x);

        for (int j = 0; j < nm_use; ++j) {
            const auto& xj = X_cols[j];
            for (const auto& ij : xj) {
                x_trips.emplace_back(ij.first, j, ij.second);
            }
        }
        X.setFromTriplets(x_trips.begin(), x_trips.end());
        X.makeCompressed();
    }

    return { std::move(T_full), std::move(X), std::move(masters_glob), std::move(slaves_glob) };
}

}} // namespace fem::constraint
