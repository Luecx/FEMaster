/**
 * @file null_space.cpp
 * @brief Implements affine null-space construction for linear constraints.
 *
 * The constraint system `C u = d` is reduced by direct substitution before a
 * sparse QR factorization constructs `u = u_p + T q`.
 *
 * @see src/constraints/transformer/null_space.h
 * @see src/constraints/transformer/constraint_map.h
 * @author Finn Eggers
 * @date 14.07.2026
 */

#include "null_space.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <vector>

namespace fem {
namespace constraint {

std::pair<ConstraintMap, ConstraintBuildReport>
build_null_space(const ConstraintSystem& system, const NullSpaceOptions& options) {
    struct PreprocessedSystem {
        SparseMatrix  C{};
        DynamicVector d{};

        std::vector<int>       used_columns{};
        std::vector<char>      fixed_columns{};
        std::vector<Precision> fixed_values{};
        std::vector<char>      kept_rows{};
    };

    ConstraintMap         map{};
    ConstraintBuildReport report{};

    report.equations   = system.equations;
    report.dofs        = system.dofs;
    report.rhs_norm    = system.d.size() == 0 ? Precision(0) : system.d.norm();
    report.homogeneous = system.d.size() == 0 ||
                         system.d.lpNorm<Eigen::Infinity>() == Precision(0);

    Time preprocess_time = 0;
    Time qr_time         = 0;
    Time rank_time       = 0;
    Time build_x_time    = 0;
    Time partition_time  = 0;
    Time assemble_time   = 0;
    Time particular_time = 0;
#ifndef NDEBUG
    Time invariants_time = 0;
#endif

    auto log_timing = [](const char* label, Time duration) {
        logging::info(true,
                      "  ", std::left, std::setw(18), label,
                      " ", std::right, std::setw(8), duration);
    };

    // Extract directly prescribed DOFs and compact unresolved columns for QR
    PreprocessedSystem preprocessed{};
    preprocess_time = Timer::measure_time([&]() {
        preprocessed.C = system.C;
        preprocessed.d = system.d.size() == system.equations
                             ? system.d
                             : DynamicVector::Zero(system.equations);

        // Identify rows that prescribe exactly one full-system DOF
        std::vector<int>       row_nonzeros(system.equations, 0);
        std::vector<int>       row_columns(system.equations, -1);
        std::vector<Precision> row_values(system.equations, Precision(0));

        for (int column = 0; column < system.C.outerSize(); ++column) {
            for (SparseMatrix::InnerIterator entry(system.C, column); entry; ++entry) {
                std::size_t row = static_cast<std::size_t>(entry.row());
                row_nonzeros[row]++;
                row_columns [row] = column;
                row_values  [row]  = entry.value();
            }
        }

        preprocessed.fixed_columns.assign(system.dofs, 0);
        preprocessed.fixed_values.assign(system.dofs, Precision(0));
        preprocessed.kept_rows.assign(system.equations, 1);

        for (Index row = 0; row < system.equations; ++row) {
            if (row_nonzeros[static_cast<std::size_t>(row)] != 1) {
                continue;
            }

            const int       column      = row_columns[static_cast<std::size_t>(row)];
            const Precision coefficient = row_values[static_cast<std::size_t>(row)];
            if (column < 0 || coefficient == Precision(0)) {
                continue;
            }

            preprocessed.fixed_columns[static_cast<std::size_t>(column)] = 1;
            preprocessed.fixed_values[static_cast<std::size_t>(column)]  =
                report.homogeneous ? Precision(0) : preprocessed.d[row] / coefficient;
            preprocessed.kept_rows[static_cast<std::size_t>(row)] = 0;
        }

        // Substitute prescribed values into all unresolved equations
        for (Index column = 0; column < system.dofs; ++column) {
            if (!preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                continue;
            }

            const Precision fixed_value = preprocessed.fixed_values[static_cast<std::size_t>(column)];
            if (fixed_value == Precision(0)) {
                continue;
            }

            for (SparseMatrix::InnerIterator entry(preprocessed.C, column); entry; ++entry) {
                preprocessed.d[entry.row()] -= entry.value() * fixed_value;
            }
        }

        // Remove prescribed columns and equations resolved by substitution
        TripletList unfixed_entries{};
        unfixed_entries.reserve(preprocessed.C.nonZeros());
        for (int column = 0; column < preprocessed.C.outerSize(); ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                continue;
            }

            for (SparseMatrix::InnerIterator entry(preprocessed.C, column); entry; ++entry) {
                if (preprocessed.kept_rows[static_cast<std::size_t>(entry.row())]) {
                    unfixed_entries.emplace_back(entry.row(), column, entry.value());
                }
            }
        }

        preprocessed.C.setZero();
        preprocessed.C.resize(system.equations, system.dofs);
        preprocessed.C.setFromTriplets(unfixed_entries.begin(), unfixed_entries.end());
        preprocessed.C.makeCompressed();

        // Compact nonzero columns while retaining their full-system indices
        for (int column = 0; column < preprocessed.C.outerSize(); ++column) {
            if (preprocessed.C.col(column).nonZeros() > 0) {
                preprocessed.used_columns.push_back(column);
            }
        }

        TripletList compact_entries{};
        compact_entries.reserve(preprocessed.C.nonZeros());
        for (int local_column = 0;
             local_column < static_cast<int>(preprocessed.used_columns.size());
             ++local_column) {
            const int global_column = preprocessed.used_columns[static_cast<std::size_t>(local_column)];
            for (SparseMatrix::InnerIterator entry(preprocessed.C, global_column); entry; ++entry) {
                compact_entries.emplace_back(entry.row(), local_column, entry.value());
            }
        }

        SparseMatrix compact(system.equations, static_cast<int>(preprocessed.used_columns.size()));
        compact.setFromTriplets(compact_entries.begin(), compact_entries.end());
        compact.makeCompressed();
        preprocessed.C = std::move(compact);
    });

    const Index n_fixed = static_cast<Index>(std::count(
        preprocessed.fixed_columns.begin(),
        preprocessed.fixed_columns.end(),
        char(1)
    ));

    // Construct the map directly when no unresolved constraints remain
    if (preprocessed.used_columns.empty()) {
        map.full_size  = system.dofs;
        map.particular = DynamicVector::Zero(system.dofs);

        for (Index column = 0; column < system.dofs; ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                map.slaves.push_back(column);
                map.particular[column] = preprocessed.fixed_values[static_cast<std::size_t>(column)];
            } else {
                map.masters.push_back(column);
            }
        }

        TripletList identity_entries{};
        identity_entries.reserve(map.masters.size());
        for (Index column = 0; column < map.n_master(); ++column) {
            identity_entries.emplace_back(map.masters[static_cast<std::size_t>(column)], column, Precision(1));
        }

        map.T.resize(system.dofs, map.n_master());
        map.T.setFromTriplets(identity_entries.begin(), identity_entries.end());
        map.T.makeCompressed();
        map.X.resize(map.n_slave(), map.n_master());

        report.rank           = n_fixed;
        report.redundant_rows = system.equations - report.rank;
        report.residual_norm  = (system.C * map.particular - system.d).norm();
        report.feasible       = std::isfinite(report.residual_norm) &&
                                report.residual_norm <= options.feasibility_tolerance *
                                    std::max<Precision>(report.rhs_norm, Precision(1));

        logging::info(true, "");
        logging::info(true, "Null-space timings (ms):");
        log_timing("Preprocess", preprocess_time);
        return {std::move(map), report};
    }

    // Factorize unresolved constraints with column pivoting
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr{};
    SparseMatrix R{};
    qr_time = Timer::measure_time([&]() {
        const Precision pivot_tolerance = std::min<Precision>(
            Precision(0.1),
            std::max<Precision>(Precision(0), options.rank_tolerance)
        );
        qr.setPivotThreshold(pivot_tolerance);
        qr.compute(preprocessed.C);
        logging::error(qr.info() == Eigen::Success, "[NullSpace] SparseQR failed");
        R = qr.matrixR();
    });

    // Determine numerical rank from the leading diagonal of R
    int       reduced_rank = 0;
    Precision max_diagonal = Precision(0);
    rank_time = Timer::measure_time([&]() {
        const int max_rank = std::min<int>(R.rows(), R.cols());
        for (int diagonal = 0; diagonal < max_rank; ++diagonal) {
            max_diagonal = std::max(max_diagonal, std::abs(R.coeff(diagonal, diagonal)));
        }

        const Precision tolerance = std::max<Precision>(
            Precision(1e-300),
            max_diagonal * options.rank_tolerance
        );
        while (reduced_rank < max_rank &&
               std::abs(R.coeff(reduced_rank, reduced_rank)) > tolerance) {
            ++reduced_rank;
        }
    });

    // Compute the slave-to-master block X = -R11^-1 R12
    using SparseColumns = std::vector<std::vector<std::pair<int, Precision>>>;
    SparseColumns x_columns{};
    build_x_time = Timer::measure_time([&]() {
        const int master_columns = R.cols() - reduced_rank;
        x_columns.resize(static_cast<std::size_t>(std::max(0, master_columns)));
        if (reduced_rank == 0 || master_columns <= 0) {
            return;
        }

        std::vector<Precision> diagonal(static_cast<std::size_t>(reduced_rank), Precision(0));
        SparseColumns          upper_rows(static_cast<std::size_t>(reduced_rank));

        // Store R11 by rows for independent sparse back substitutions
        for (int column = 0; column < reduced_rank; ++column) {
            for (SparseMatrix::InnerIterator entry(R, column); entry; ++entry) {
                const int row = entry.row();
                if (row >= reduced_rank) {
                    continue;
                }
                if (row == column) {
                    diagonal[static_cast<std::size_t>(row)] = entry.value();
                } else if (row < column) {
                    upper_rows[static_cast<std::size_t>(row)].emplace_back(column, entry.value());
                }
            }
        }

        // Each master column defines an independent triangular solve
#ifdef _OPENMP
#pragma omp parallel for if (master_columns > 2)
#endif
        for (int master = 0; master < master_columns; ++master) {
            std::unordered_map<int, Precision> rhs{};
            rhs.reserve(16);
            for (SparseMatrix::InnerIterator entry(R, reduced_rank + master); entry; ++entry) {
                if (entry.row() < reduced_rank) {
                    rhs[entry.row()] = -entry.value();
                }
            }
            if (rhs.empty()) {
                continue;
            }

            std::vector<Precision> values(static_cast<std::size_t>(reduced_rank), Precision(0));
            for (int row = reduced_rank - 1; row >= 0; --row) {
                Precision value = Precision(0);
                const auto found = rhs.find(row);
                if (found != rhs.end()) {
                    value = found->second;
                }
                for (const auto& [column, coefficient] : upper_rows[static_cast<std::size_t>(row)]) {
                    value -= coefficient * values[static_cast<std::size_t>(column)];
                }

                const Precision diagonal_value = diagonal[static_cast<std::size_t>(row)];
                values[static_cast<std::size_t>(row)] =
                    diagonal_value != Precision(0) ? value / diagonal_value : Precision(0);
            }

            auto& entries = x_columns[static_cast<std::size_t>(master)];
            entries.reserve(16);
            for (int row = 0; row < reduced_rank; ++row) {
                const Precision value = values[static_cast<std::size_t>(row)];
                if (value != Precision(0)) {
                    entries.emplace_back(row, value);
                }
            }
        }
    });

    // Map the QR partition back to full-system DOF indices
    std::vector<int> qr_slaves{};
    std::vector<int> qr_masters{};
    partition_time = Timer::measure_time([&]() {
        const auto& permutation = qr.colsPermutation().indices();
        qr_slaves.reserve(static_cast<std::size_t>(reduced_rank));
        for (int row = 0; row < reduced_rank; ++row) {
            qr_slaves.push_back(permutation(row));
        }

        const int reduced_masters = static_cast<int>(preprocessed.used_columns.size()) - reduced_rank;
        qr_masters.reserve(static_cast<std::size_t>(std::max(0, reduced_masters)));
        for (int column = 0; column < reduced_masters; ++column) {
            qr_masters.push_back(permutation(reduced_rank + column));
        }

        std::vector<char> used_global(static_cast<std::size_t>(system.dofs), 0);
        for (int column : preprocessed.used_columns) {
            used_global[static_cast<std::size_t>(column)] = 1;
        }

        map.full_size = system.dofs;
        map.slaves.reserve(static_cast<std::size_t>(n_fixed + reduced_rank));
        for (Index column = 0; column < system.dofs; ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                map.slaves.push_back(column);
            }
        }
        for (int local_column : qr_slaves) {
            map.slaves.push_back(preprocessed.used_columns[static_cast<std::size_t>(local_column)]);
        }

        map.masters.reserve(static_cast<std::size_t>(system.dofs - n_fixed - reduced_rank));
        for (int local_column : qr_masters) {
            map.masters.push_back(preprocessed.used_columns[static_cast<std::size_t>(local_column)]);
        }
        for (Index column = 0; column < system.dofs; ++column) {
            if (!used_global[static_cast<std::size_t>(column)] &&
                !preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                map.masters.push_back(column);
            }
        }
    });

    // Assemble the affine map u = u_p + T q
    assemble_time = Timer::measure_time([&]() {
        map.particular = DynamicVector::Zero(system.dofs);
        for (Index column = 0; column < system.dofs; ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                map.particular[column] = preprocessed.fixed_values[static_cast<std::size_t>(column)];
            }
        }

        std::unordered_map<Index, Index> master_local{};
        master_local.reserve(map.masters.size());
        for (Index column = 0; column < map.n_master(); ++column) {
            master_local[map.masters[static_cast<std::size_t>(column)]] = column;
        }

        TripletList t_entries{};
        TripletList x_entries{};
        t_entries.reserve(map.masters.size() + R.nonZeros());
        x_entries.reserve(R.nonZeros());

        for (Index column = 0; column < map.n_master(); ++column) {
            t_entries.emplace_back(map.masters[static_cast<std::size_t>(column)], column, Precision(1));
        }

        for (int master = 0; master < static_cast<int>(qr_masters.size()); ++master) {
            const Index global_master = preprocessed.used_columns[
                static_cast<std::size_t>(qr_masters[static_cast<std::size_t>(master)])
            ];
            const Index map_column = master_local.at(global_master);

            for (const auto& [slave_row, coefficient] : x_columns[static_cast<std::size_t>(master)]) {
                const Index global_slave = preprocessed.used_columns[
                    static_cast<std::size_t>(qr_slaves[static_cast<std::size_t>(slave_row)])
                ];
                const Index map_row = n_fixed + slave_row;
                t_entries.emplace_back(global_slave, map_column, coefficient);
                x_entries.emplace_back(map_row, map_column, coefficient);
            }
        }

        map.T.resize(system.dofs, map.n_master());
        map.T.setFromTriplets(t_entries.begin(), t_entries.end());
        map.T.makeCompressed();

        map.X.resize(map.n_slave(), map.n_master());
        map.X.setFromTriplets(x_entries.begin(), x_entries.end());
        map.X.makeCompressed();
    });

    // Compute a particular solution for inhomogeneous constraints
    particular_time = Timer::measure_time([&]() {
        if (report.homogeneous) {
            report.residual_norm = (system.C * map.particular - system.d).norm();
            report.feasible      = true;
            return;
        }

        const DynamicVector reduced_solution = qr.solve(preprocessed.d);
        logging::error(qr.info() == Eigen::Success, "[NullSpace] particular solution failed");

        DynamicVector full_solution = DynamicVector::Zero(system.dofs);
        for (int local_column = 0;
             local_column < static_cast<int>(preprocessed.used_columns.size());
             ++local_column) {
            full_solution[preprocessed.used_columns[static_cast<std::size_t>(local_column)]] =
                reduced_solution[local_column];
        }
        for (Index column = 0; column < system.dofs; ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                full_solution[column] = preprocessed.fixed_values[static_cast<std::size_t>(column)];
            }
        }

        DynamicVector master_values(map.n_master());
        for (Index column = 0; column < map.n_master(); ++column) {
            master_values[column] = full_solution[map.masters[static_cast<std::size_t>(column)]];
        }

        const DynamicVector slave_values = map.X * master_values;

        // Remove the homogeneous master contribution from the slave values
        for (Index row = 0; row < map.n_slave(); ++row) {
            const Index global_slave = map.slaves[static_cast<std::size_t>(row)];
            map.particular[global_slave] = full_solution[global_slave] - slave_values[row];
        }
        for (Index column = 0; column < system.dofs; ++column) {
            if (preprocessed.fixed_columns[static_cast<std::size_t>(column)]) {
                map.particular[column] = preprocessed.fixed_values[static_cast<std::size_t>(column)];
            }
        }

        report.residual_norm = (system.C * map.particular - system.d).norm();
        report.feasible      = report.residual_norm <= options.feasibility_tolerance *
                               std::max<Precision>(report.rhs_norm, Precision(1));
    });

    report.rank           = n_fixed + reduced_rank;
    report.redundant_rows = std::max<Index>(0, system.equations - report.rank);

#ifndef NDEBUG
    // Verify that every transformation column lies in the null space of C
    invariants_time = Timer::measure_time([&]() {
        Precision max_null_space_norm = Precision(0);
        for (Index column = 0; column < map.n_master(); ++column) {
            max_null_space_norm = std::max(
                max_null_space_norm,
                (system.C * map.T.col(column)).norm()
            );
        }
        logging::error(max_null_space_norm <= Precision(1e-12),
                       "[NullSpace] invariant failed: max ||C*T(:,j)|| = ",
                       max_null_space_norm);
    });
#endif

    logging::info(true, "");
    logging::info(true, "Null-space timings (ms):");
    log_timing("Preprocess", preprocess_time);
    log_timing("QR",         qr_time);
    log_timing("Rank",       rank_time);
    log_timing("Build X",    build_x_time);
    log_timing("Partition",  partition_time);
    log_timing("Assemble",   assemble_time);
#ifndef NDEBUG
    log_timing("Invariants", invariants_time);
#endif
    log_timing("Particular", particular_time);

    return {std::move(map), report};
}

} // namespace constraint
} // namespace fem
