/**
 * @file elimination.cpp
 * @brief Implements direct substitution of linear constraint equations.
 *
 * Each independent equation selects one dependent DOF. Existing substitutions
 * are applied before pivot selection to keep the dependency graph acyclic.
 *
 * @see src/constraints/transformer/elimination.h
 * @see src/constraints/transformer/constraint_map.h
 * @author Finn Eggers
 * @date 14.07.2026
 */

#include "elimination.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>

namespace fem {
namespace constraint {

std::pair<ConstraintMap, ConstraintBuildReport> build_elimination(const ConstraintSystem& system) {
    constexpr Precision tolerance = Precision(1e-14);

    ConstraintMap         map{};
    ConstraintBuildReport report{};

    report.equations   = system.equations;
    report.dofs        = system.dofs;
    report.homogeneous = system.d.size() == 0 ||
                         system.d.lpNorm<Eigen::Infinity>() == Precision(0);
    report.rhs_norm    = system.d.size() == 0 ? Precision(0) : system.d.norm();

    // Store equations by row and count each DOF's global participation
    std::vector<std::vector<std::pair<Index, Precision>>> rows(
        static_cast<std::size_t>(system.equations)
    );
    std::vector<Index> column_count(static_cast<std::size_t>(system.dofs), 0);

    for (int column = 0; column < system.C.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(system.C, column); entry; ++entry) {
            if (std::abs(entry.value()) <= tolerance) {
                continue;
            }

            rows[static_cast<std::size_t>(entry.row())].emplace_back(column, entry.value());
            ++column_count[static_cast<std::size_t>(column)];
        }
    }

    std::vector<char> is_slave(static_cast<std::size_t>(system.dofs), 0);
    std::vector<Precision> pivot_constant(static_cast<std::size_t>(system.dofs), Precision(0));
    std::vector<std::vector<std::pair<Index, Precision>>> pivot_dependencies(
        static_cast<std::size_t>(system.dofs)
    );

    map.slaves.reserve(static_cast<std::size_t>(system.equations));

    struct Expanded {
        Precision constant = Precision(0);
        std::unordered_map<Index, Precision> coefficients{};
    };

    // Express a DOF through the master DOFs available at the current row
    auto expand_current = [&](auto&& self, Index column) -> Expanded {
        if (!is_slave[static_cast<std::size_t>(column)]) {
            Expanded result{};
            result.coefficients[column] = Precision(1);
            return result;
        }

        Expanded result{};
        result.constant = pivot_constant[static_cast<std::size_t>(column)];
        for (const auto& [dependency, coefficient] :
             pivot_dependencies[static_cast<std::size_t>(column)]) {
            const Expanded expanded_dependency = self(self, dependency);
            result.constant += coefficient * expanded_dependency.constant;
            for (const auto& [master, value] : expanded_dependency.coefficients) {
                result.coefficients[master] += coefficient * value;
            }
        }
        return result;
    };

    // Eliminate previously selected slave DOFs before choosing each new pivot
    for (Index row = 0; row < system.equations; ++row) {
        const auto& entries = rows[static_cast<std::size_t>(row)];
        Precision rhs = system.d.size() == system.equations ? system.d[row] : Precision(0);

        if (entries.empty()) {
            logging::error(std::abs(rhs) <= tolerance,
                           "[Elimination] inconsistent empty constraint row");
            continue;
        }

        std::unordered_map<Index, Precision> reduced_coefficients{};
        reduced_coefficients.reserve(entries.size() * 2);
        for (const auto& [column, coefficient] : entries) {
            const Expanded expanded_column = expand_current(expand_current, column);
            rhs -= coefficient * expanded_column.constant;
            for (const auto& [expanded_column_id, expanded_value] : expanded_column.coefficients) {
                reduced_coefficients[expanded_column_id] += coefficient * expanded_value;
            }
        }

        std::vector<std::pair<Index, Precision>> reduced_entries{};
        reduced_entries.reserve(reduced_coefficients.size());
        for (const auto& [column, coefficient] : reduced_coefficients) {
            if (std::abs(coefficient) > tolerance) {
                reduced_entries.emplace_back(column, coefficient);
            }
        }

        if (reduced_entries.empty()) {
            logging::error(std::abs(rhs) <= tolerance,
                           "[Elimination] inconsistent redundant constraint row");
            ++report.redundant_rows;
            continue;
        }

        // Prefer sparsely occurring coefficients close to unit magnitude
        Index     pivot_column      = system.dofs;
        Precision pivot_coefficient = Precision(0);
        Index     best_count        = std::numeric_limits<Index>::max();
        Precision best_unit_distance = std::numeric_limits<Precision>::max();

        for (const auto& [column, coefficient] : reduced_entries) {
            if (is_slave[static_cast<std::size_t>(column)]) {
                continue;
            }

            const Index     count         = column_count[static_cast<std::size_t>(column)];
            const Precision unit_distance = std::abs(std::abs(coefficient) - Precision(1));
            const bool better =
                pivot_column == system.dofs ||
                count < best_count ||
                (count == best_count && unit_distance < best_unit_distance) ||
                (count == best_count &&
                 unit_distance == best_unit_distance &&
                 std::abs(coefficient) > std::abs(pivot_coefficient));

            if (better) {
                pivot_column       = column;
                pivot_coefficient  = coefficient;
                best_count         = count;
                best_unit_distance = unit_distance;
            }
        }

        if (pivot_column == system.dofs) {
            ++report.redundant_rows;
            continue;
        }

        logging::error(std::abs(pivot_coefficient) > tolerance,
                       "[Elimination] pivot is numerically zero");

        is_slave[static_cast<std::size_t>(pivot_column)] = 1;
        map.slaves.push_back(pivot_column);
        pivot_constant[static_cast<std::size_t>(pivot_column)] = rhs / pivot_coefficient;

        // Store u_pivot as an affine combination of the remaining DOFs
        auto& dependencies = pivot_dependencies[static_cast<std::size_t>(pivot_column)];
        dependencies.reserve(reduced_entries.size());
        for (const auto& [column, coefficient] : reduced_entries) {
            if (column != pivot_column) {
                dependencies.emplace_back(column, -coefficient / pivot_coefficient);
            }
        }
    }

    // Collect independent DOFs in full-system order
    map.full_size = system.dofs;
    map.masters.reserve(system.dofs - map.n_slave());
    for (Index column = 0; column < system.dofs; ++column) {
        if (!is_slave[column]) {
            map.masters.push_back(column);
        }
    }

    std::unordered_map<Index, Index> master_local{};
    master_local.reserve(map.masters.size() * 2);
    for (Index column = 0; column < map.n_master(); ++column) {
        master_local[map.masters[column]] = column;
    }

    std::vector<char>     expansion_state(system.dofs, 0);
    std::vector<Expanded> expanded_values(system.dofs);

    // Recursively express every slave through final master DOFs
    std::function<Expanded(Index)> expand = [&](Index column) -> Expanded {
        if (!is_slave[column]) {
            Expanded result{};
            result.coefficients[column] = Precision(1);
            return result;
        }

        char& state = expansion_state[column];
        logging::error(state != 1, "[Elimination] cyclic constraint dependency");
        if (state == 2) {
            return expanded_values[column];
        }

        state = 1;
        Expanded result{};
        result.constant = pivot_constant[column];
        for (const auto& [dependency, coefficient] :
             pivot_dependencies[column]) {
            const Expanded expanded_dependency = expand(dependency);
            result.constant += coefficient * expanded_dependency.constant;
            for (const auto& [master, value] : expanded_dependency.coefficients) {
                result.coefficients[master] += coefficient * value;
            }
        }

        state = 2;
        expanded_values[static_cast<std::size_t>(column)] = result;
        return result;
    };

    // Assemble u = u_p + T q and the explicit slave block X
    TripletList t_entries{};
    TripletList x_entries{};
    t_entries.reserve(map.masters.size() + system.C.nonZeros());
    x_entries.reserve(system.C.nonZeros());

    map.particular = DynamicVector::Zero(system.dofs);
    for (Index column = 0; column < map.n_master(); ++column) {
        t_entries.emplace_back(map.masters[static_cast<std::size_t>(column)], column, Precision(1));
    }

    for (Index row = 0; row < map.n_slave(); ++row) {
        const Index    slave = map.slaves[static_cast<std::size_t>(row)];
        const Expanded value = expand(slave);
        map.particular[slave] = value.constant;

        for (const auto& [master, coefficient] : value.coefficients) {
            if (std::abs(coefficient) <= tolerance) {
                continue;
            }

            const auto found = master_local.find(master);
            logging::error(found != master_local.end(),
                           "[Elimination] dependency expanded to a non-master DOF");

            t_entries.emplace_back(slave, found->second, coefficient);
            x_entries.emplace_back(row, found->second, coefficient);
        }
    }

    map.T.resize(system.dofs, map.n_master());
    map.T.setFromTriplets(t_entries.begin(), t_entries.end());
    map.T.makeCompressed();

    map.X.resize(map.n_slave(), map.n_master());
    map.X.setFromTriplets(x_entries.begin(), x_entries.end());
    map.X.makeCompressed();

    // Evaluate feasibility against the original assembled constraints
    report.rank          = map.n_slave();
    report.residual_norm = (system.C * map.particular - system.d).norm();
    report.feasible      = std::isfinite(report.residual_norm) &&
                           report.residual_norm <= Precision(1e-10) *
                               std::max<Precision>(report.rhs_norm, Precision(1));

    logging::error(std::isfinite(report.residual_norm),
                   "[Elimination] invalid particular solution");
    return {std::move(map), report};
}

} // namespace constraint
} // namespace fem
