#include "amg_preconditioner.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace fem::solver::detail {

namespace {
constexpr Precision diagonal_epsilon = 1e-30;
constexpr Precision spectral_safety = 1.25;

Precision coupling_strength(Precision value,
                            Precision inv_diag_i,
                            Precision inv_diag_j) {
    return std::abs(value) * std::sqrt(std::abs(inv_diag_i * inv_diag_j));
}
} // namespace

void CpuAmgPreconditioner::compute(const SparseMatrix& matrix) {
    logging::error(matrix.rows() == matrix.cols(), "AMG requires a square matrix");
    logging::error(matrix.rows() > 0, "AMG requires a non-empty matrix");

    logging::info(true, "Building CPU local-spectral aggregation AMG hierarchy");
    logging::up();
    logging::info(true, "fine level: N=", matrix.rows(), " nnz=", matrix.nonZeros());

    Timer timer {};
    timer.start();

    SparseMatrix scaled = matrix;
    m_scaling = equilibrate(scaled);
    scaled.makeCompressed();

    timer.stop();
    logging::info(true, "symmetric diagonal scaling: ", timer.elapsed(), " ms");

    m_levels.clear();
    m_workspace.clear();

    Level first;
    first.a = std::move(scaled);
    m_levels.push_back(std::move(first));

    while (static_cast<int>(m_levels.size()) < max_levels
           && m_levels.back().a.rows() > coarse_direct_size) {
        const std::size_t level_index = m_levels.size() - 1;
        Level& fine = m_levels.back();

        timer.start();
        prepare_level(fine);
        timer.stop();
        logging::info(true,
                      "level ", level_index,
                      " smoother setup: ", timer.elapsed(), " ms",
                      " omega=", fine.smoother_omega);

        timer.start();
        int aggregate_count = 0;
        int min_aggregate_size = 0;
        int max_aggregate_size = 0;
        Precision average_aggregate_size = 0;
        const auto aggregates = build_aggregates(
            fine.a,
            aggregate_count,
            min_aggregate_size,
            max_aggregate_size,
            average_aggregate_size);
        timer.stop();

        logging::info(true,
                      "level ", level_index,
                      " algebraic aggregation: ", timer.elapsed(), " ms",
                      " -> aggregates=", aggregate_count,
                      " size[min/avg/max]=", min_aggregate_size,
                      "/", average_aggregate_size,
                      "/", max_aggregate_size);

        if (aggregate_count <= 0 || aggregate_count >= fine.a.rows()) {
            logging::warning(false, "AMG aggregation did not reduce the system; stopping hierarchy construction");
            break;
        }

        timer.start();
        Eigen::Index coarse_size = 0;
        std::vector<Eigen::Index> rank_histogram(static_cast<std::size_t>(max_local_basis + 1), 0);
        SparseMatrix prolongator = build_spectral_prolongator(
            fine.a,
            aggregates,
            aggregate_count,
            coarse_size,
            rank_histogram);
        timer.stop();

        logging::info(true,
                      "level ", level_index,
                      " local spectral basis: ", timer.elapsed(), " ms",
                      " -> N=", coarse_size,
                      " P nnz=", prolongator.nonZeros());
        logging::info(true,
                      "level ", level_index,
                      " local ranks: 1=", rank_histogram[1],
                      " 2=", rank_histogram[2],
                      " 3=", rank_histogram[3],
                      " 4=", rank_histogram[4],
                      " 5=", rank_histogram[5],
                      " 6=", rank_histogram[6]);

        if (coarse_size <= 0 || coarse_size >= fine.a.rows()) {
            logging::warning(false, "AMG spectral coarse space did not reduce the system; stopping hierarchy construction");
            break;
        }

        timer.start();
        SparseMatrix ap = fine.a * prolongator;
        ap.makeCompressed();
        timer.stop();
        logging::info(true,
                      "level ", level_index,
                      " A*P: ", timer.elapsed(), " ms",
                      " -> nnz=", ap.nonZeros());

        timer.start();
        SparseMatrix prolongator_t = prolongator.transpose();
        prolongator_t.makeCompressed();
        SparseMatrix coarse = prolongator_t * ap;
        coarse.makeCompressed();
        timer.stop();
        logging::info(true,
                      "level ", level_index,
                      " P^T*A*P: ", timer.elapsed(), " ms",
                      " -> N=", coarse.rows(),
                      " nnz=", coarse.nonZeros());

        timer.start();
        const DynamicVector coarse_scaling = equilibrate(coarse);
        scale_prolongator_columns(prolongator, coarse_scaling);
        prolongator.makeCompressed();
        prolongator_t = prolongator.transpose();
        prolongator_t.makeCompressed();
        timer.stop();
        logging::info(true,
                      "level ", level_index,
                      " coarse equilibration: ", timer.elapsed(), " ms");

        fine.p = std::move(prolongator);
        fine.pt = std::move(prolongator_t);

        Level next;
        next.a = std::move(coarse);
        m_levels.push_back(std::move(next));
    }

    timer.start();
    prepare_level(m_levels.back());
    timer.stop();
    logging::info(true,
                  "coarsest level: N=", m_levels.back().a.rows(),
                  " nnz=", m_levels.back().a.nonZeros(),
                  " smoother setup=", timer.elapsed(), " ms",
                  " omega=", m_levels.back().smoother_omega);

    m_workspace.resize(m_levels.size());
    for (std::size_t level = 0; level + 1 < m_levels.size(); ++level) {
        m_workspace[level].residual.resize(m_levels[level].a.rows());
        m_workspace[level].coarse_rhs.resize(m_levels[level + 1].a.rows());
        m_workspace[level].coarse_error.resize(m_levels[level + 1].a.rows());
    }

    m_scaled_rhs.resize(matrix.rows());
    m_scaled_solution.resize(matrix.rows());

    timer.start();
    m_coarse_solver.compute(m_levels.back().a);
    timer.stop();
    logging::error(m_coarse_solver.info() == Eigen::Success,
                   "AMG coarse-level LDLT factorization failed");
    logging::info(true, "coarse LDLT: ", timer.elapsed(), " ms");
    logging::down();
}

void CpuAmgPreconditioner::apply(const DynamicVector& rhs,
                                 DynamicVector& solution) const {
    logging::error(rhs.size() == m_scaling.size(), "AMG RHS size mismatch");

    m_scaled_rhs.array() = m_scaling.array() * rhs.array();
    m_scaled_solution.setZero();
    vcycle(0, m_scaled_rhs, m_scaled_solution);

    solution.resize(rhs.size());
    solution.array() = m_scaling.array() * m_scaled_solution.array();
}

DynamicVector CpuAmgPreconditioner::inverse_diagonal(const SparseMatrix& matrix) {
    DynamicVector diagonal = DynamicVector::Zero(matrix.rows());

    for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry) {
            if (entry.row() == entry.col()) {
                diagonal[entry.row()] = entry.value();
            }
        }
    }

    DynamicVector inv_diag(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        logging::error(diagonal[i] > diagonal_epsilon,
                       "AMG requires a positive diagonal, invalid entry at row ", i, ": ", diagonal[i]);
        inv_diag[i] = Precision(1) / diagonal[i];
    }
    return inv_diag;
}

DynamicVector CpuAmgPreconditioner::equilibrate(SparseMatrix& matrix) {
    const DynamicVector inv_diag = inverse_diagonal(matrix);
    DynamicVector scaling(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        scaling[i] = std::sqrt(inv_diag[i]);
    }

    for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry) {
            entry.valueRef() *= scaling[entry.row()] * scaling[entry.col()];
        }
    }

    return scaling;
}

Precision CpuAmgPreconditioner::estimate_spectral_radius(
        const SparseMatrix& matrix,
        const DynamicVector& inv_diag) {
    DynamicVector inv_sqrt(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        inv_sqrt[i] = std::sqrt(inv_diag[i]);
    }

    DynamicVector x(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        const Precision index = static_cast<Precision>(i + 1);
        x[i] = std::sin(index * Precision(12.9898))
             + Precision(0.5) * std::sin(index * Precision(78.233));
    }

    const Precision initial_norm = x.norm();
    if (!(initial_norm > diagonal_epsilon) || !std::isfinite(initial_norm)) {
        return Precision(1);
    }
    x /= initial_norm;

    Precision rho = Precision(1);
    for (int iteration = 0; iteration < spectral_iterations; ++iteration) {
        const DynamicVector scaled_x = inv_sqrt.array() * x.array();
        DynamicVector y = matrix * scaled_x;
        y.array() *= inv_sqrt.array();

        const Precision norm = y.norm();
        if (!(norm > diagonal_epsilon) || !std::isfinite(norm)) {
            return Precision(1);
        }

        x = y / norm;
        rho = norm;
    }

    return std::max(Precision(1), spectral_safety * rho);
}

std::vector<int> CpuAmgPreconditioner::pair_graph(
        const SparseMatrix& graph,
        const std::vector<int>& group_sizes,
        int& coarse_count,
        std::vector<int>& coarse_sizes) {
    const int n = static_cast<int>(graph.rows());
    const DynamicVector inv_diag = inverse_diagonal(graph);

    std::vector<int> aggregates(static_cast<std::size_t>(n), -1);
    coarse_sizes.clear();
    coarse_sizes.reserve(static_cast<std::size_t>((n + 1) / 2));
    coarse_count = 0;

    for (int i = 0; i < n; ++i) {
        if (aggregates[static_cast<std::size_t>(i)] >= 0) {
            continue;
        }

        int strongest_neighbor = -1;
        Precision strongest_strength = -std::numeric_limits<Precision>::infinity();

        for (SparseMatrix::InnerIterator entry(graph, i); entry; ++entry) {
            const int j = static_cast<int>(entry.row());
            if (j == i || aggregates[static_cast<std::size_t>(j)] >= 0) {
                continue;
            }
            if (group_sizes[static_cast<std::size_t>(i)]
                + group_sizes[static_cast<std::size_t>(j)] > target_aggregate_size) {
                continue;
            }

            const Precision strength = coupling_strength(entry.value(), inv_diag[i], inv_diag[j]);
            if (strength > strongest_strength) {
                strongest_strength = strength;
                strongest_neighbor = j;
            }
        }

        const int aggregate = coarse_count++;
        aggregates[static_cast<std::size_t>(i)] = aggregate;
        int size = group_sizes[static_cast<std::size_t>(i)];

        if (strongest_neighbor >= 0 && strongest_strength > Precision(0)) {
            aggregates[static_cast<std::size_t>(strongest_neighbor)] = aggregate;
            size += group_sizes[static_cast<std::size_t>(strongest_neighbor)];
        }

        coarse_sizes.push_back(size);
    }

    return aggregates;
}

SparseMatrix CpuAmgPreconditioner::aggregate_strength_graph(
        const SparseMatrix& graph,
        const std::vector<int>& aggregates,
        int aggregate_count) {
    const int n = static_cast<int>(graph.rows());

    std::vector<int> first(static_cast<std::size_t>(aggregate_count), -1);
    std::vector<int> next(static_cast<std::size_t>(n), -1);
    for (int fine_column = 0; fine_column < n; ++fine_column) {
        const int coarse_column = aggregates[static_cast<std::size_t>(fine_column)];
        next[static_cast<std::size_t>(fine_column)] = first[static_cast<std::size_t>(coarse_column)];
        first[static_cast<std::size_t>(coarse_column)] = fine_column;
    }

    std::vector<int> marker(static_cast<std::size_t>(aggregate_count), -1);
    std::vector<Precision> values(static_cast<std::size_t>(aggregate_count), Precision(0));
    std::vector<int> touched;
    touched.reserve(256);

    TripletList triplets;
    const Eigen::Index reserve_count = std::min<Eigen::Index>(
        graph.nonZeros(), static_cast<Eigen::Index>(aggregate_count) * 128);
    triplets.reserve(static_cast<std::size_t>(reserve_count));

    for (int coarse_column = 0; coarse_column < aggregate_count; ++coarse_column) {
        touched.clear();

        for (int fine_column = first[static_cast<std::size_t>(coarse_column)];
             fine_column >= 0;
             fine_column = next[static_cast<std::size_t>(fine_column)]) {
            for (SparseMatrix::InnerIterator entry(graph, fine_column); entry; ++entry) {
                const int coarse_row = aggregates[static_cast<std::size_t>(entry.row())];
                const Precision value = std::abs(entry.value());

                if (marker[static_cast<std::size_t>(coarse_row)] != coarse_column) {
                    marker[static_cast<std::size_t>(coarse_row)] = coarse_column;
                    values[static_cast<std::size_t>(coarse_row)] = value;
                    touched.push_back(coarse_row);
                } else {
                    values[static_cast<std::size_t>(coarse_row)] += value;
                }
            }
        }

        std::sort(touched.begin(), touched.end());
        for (const int coarse_row : touched) {
            const Precision value = values[static_cast<std::size_t>(coarse_row)];
            if (value > Precision(0)) {
                triplets.emplace_back(coarse_row, coarse_column, value);
            }
        }
    }

    SparseMatrix coarse(aggregate_count, aggregate_count);
    coarse.setFromTriplets(triplets.begin(), triplets.end());
    coarse.makeCompressed();
    return coarse;
}

std::vector<int> CpuAmgPreconditioner::build_aggregates(
        const SparseMatrix& matrix,
        int& aggregate_count,
        int& min_aggregate_size,
        int& max_aggregate_size,
        Precision& average_aggregate_size) {
    const int n = static_cast<int>(matrix.rows());

    std::vector<int> fine_to_group(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        fine_to_group[static_cast<std::size_t>(i)] = i;
    }

    SparseMatrix graph = matrix;
    std::vector<int> group_sizes(static_cast<std::size_t>(n), 1);
    int current_count = n;

    for (int pass = 0; pass < max_aggregation_passes; ++pass) {
        int next_count = 0;
        std::vector<int> next_sizes;
        const auto pair_map = pair_graph(graph, group_sizes, next_count, next_sizes);

        if (next_count <= 0 || next_count >= current_count) {
            break;
        }

        for (int i = 0; i < n; ++i) {
            fine_to_group[static_cast<std::size_t>(i)] =
                pair_map[static_cast<std::size_t>(fine_to_group[static_cast<std::size_t>(i)])];
        }

        const int reduction = current_count - next_count;
        current_count = next_count;
        group_sizes = std::move(next_sizes);

        if (pass + 1 < max_aggregation_passes && reduction > 0) {
            graph = aggregate_strength_graph(graph, pair_map, current_count);
        }

        if (reduction <= std::max(1, current_count / 100)) {
            break;
        }
    }

    aggregate_count = current_count;
    std::vector<int> counts(static_cast<std::size_t>(aggregate_count), 0);
    for (const int aggregate : fine_to_group) {
        ++counts[static_cast<std::size_t>(aggregate)];
    }

    min_aggregate_size = n;
    max_aggregate_size = 0;
    for (const int count : counts) {
        min_aggregate_size = std::min(min_aggregate_size, count);
        max_aggregate_size = std::max(max_aggregate_size, count);
    }
    average_aggregate_size = aggregate_count > 0
        ? static_cast<Precision>(n) / static_cast<Precision>(aggregate_count)
        : Precision(0);

    return fine_to_group;
}

SparseMatrix CpuAmgPreconditioner::build_spectral_prolongator(
        const SparseMatrix& matrix,
        const std::vector<int>& aggregates,
        int aggregate_count,
        Eigen::Index& coarse_size,
        std::vector<Eigen::Index>& rank_histogram) {
    const int n = static_cast<int>(matrix.rows());

    std::vector<std::vector<int>> members(static_cast<std::size_t>(aggregate_count));
    for (int i = 0; i < n; ++i) {
        members[static_cast<std::size_t>(aggregates[static_cast<std::size_t>(i)])].push_back(i);
    }

    std::vector<int> local_position(static_cast<std::size_t>(n), -1);
    TripletList triplets;
    triplets.reserve(static_cast<std::size_t>(matrix.rows()) * static_cast<std::size_t>(max_local_basis));

    coarse_size = 0;

    for (int aggregate = 0; aggregate < aggregate_count; ++aggregate) {
        const auto& aggregate_members = members[static_cast<std::size_t>(aggregate)];
        const int local_size = static_cast<int>(aggregate_members.size());
        if (local_size == 0) {
            continue;
        }

        for (int local = 0; local < local_size; ++local) {
            local_position[static_cast<std::size_t>(aggregate_members[static_cast<std::size_t>(local)])] = local;
        }

        DynamicMatrix local = DynamicMatrix::Zero(local_size, local_size);
        for (int local_column = 0; local_column < local_size; ++local_column) {
            const int fine_column = aggregate_members[static_cast<std::size_t>(local_column)];
            for (SparseMatrix::InnerIterator entry(matrix, fine_column); entry; ++entry) {
                const int fine_row = static_cast<int>(entry.row());
                if (aggregates[static_cast<std::size_t>(fine_row)] != aggregate) {
                    continue;
                }
                const int local_row = local_position[static_cast<std::size_t>(fine_row)];
                if (local_row >= 0) {
                    local(local_row, local_column) = entry.value();
                }
            }
        }

        local = Precision(0.5) * (local + local.transpose()).eval();

        Eigen::SelfAdjointEigenSolver<DynamicMatrix> eig(local);
        logging::error(eig.info() == Eigen::Success,
                       "AMG local spectral decomposition failed for aggregate ", aggregate);

        const auto& eigenvalues = eig.eigenvalues();
        const Precision lambda_max = std::max(eigenvalues[local_size - 1], diagonal_epsilon);
        const Precision cutoff = local_spectral_cutoff * lambda_max;

        const int rank_limit = std::min(
            max_local_basis,
            std::max(1, local_size / 2));

        int rank = 1;
        for (int mode = 0; mode < rank_limit; ++mode) {
            if (eigenvalues[mode] <= cutoff) {
                rank = mode + 1;
            }
        }
        rank = std::min(rank, rank_limit);

        ++rank_histogram[static_cast<std::size_t>(rank)];

        const Eigen::Index base = coarse_size;
        for (int local_row = 0; local_row < local_size; ++local_row) {
            const int fine_row = aggregate_members[static_cast<std::size_t>(local_row)];
            for (int mode = 0; mode < rank; ++mode) {
                const Precision value = eig.eigenvectors()(local_row, mode);
                if (std::abs(value) > Precision(1e-14)) {
                    triplets.emplace_back(fine_row, base + mode, value);
                }
            }
        }

        coarse_size += rank;

        for (const int fine : aggregate_members) {
            local_position[static_cast<std::size_t>(fine)] = -1;
        }
    }

    SparseMatrix prolongator(matrix.rows(), coarse_size);
    prolongator.setFromTriplets(triplets.begin(), triplets.end());
    prolongator.makeCompressed();
    return prolongator;
}

void CpuAmgPreconditioner::scale_prolongator_columns(
        SparseMatrix& prolongator,
        const DynamicVector& scaling) {
    logging::error(prolongator.cols() == scaling.size(), "AMG prolongator scaling size mismatch");

    for (Eigen::Index column = 0; column < prolongator.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(prolongator, column); entry; ++entry) {
            entry.valueRef() *= scaling[column];
        }
    }
}

void CpuAmgPreconditioner::prepare_level(Level& level) {
    level.inv_diag = inverse_diagonal(level.a);
    const Precision spectral_radius = estimate_spectral_radius(level.a, level.inv_diag);
    level.smoother_omega = Precision(0.8) / spectral_radius;
}

void CpuAmgPreconditioner::smooth(const Level& level,
                                  const DynamicVector& rhs,
                                  DynamicVector& solution,
                                  int sweeps) const {
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        const DynamicVector residual = rhs - level.a * solution;
        solution.array() += level.smoother_omega * level.inv_diag.array() * residual.array();
    }
}

void CpuAmgPreconditioner::vcycle(std::size_t level_index,
                                  const DynamicVector& rhs,
                                  DynamicVector& solution) const {
    if (level_index + 1 == m_levels.size()) {
        solution = m_coarse_solver.solve(rhs);
        logging::error(m_coarse_solver.info() == Eigen::Success,
                       "AMG coarse-level solve failed");
        return;
    }

    const Level& level = m_levels[level_index];
    Workspace& workspace = m_workspace[level_index];

    smooth(level, rhs, solution, pre_sweeps);

    workspace.residual.noalias() = rhs - level.a * solution;
    workspace.coarse_rhs.noalias() = level.pt * workspace.residual;
    workspace.coarse_error.setZero();

    vcycle(level_index + 1, workspace.coarse_rhs, workspace.coarse_error);

    solution.noalias() += level.p * workspace.coarse_error;
    smooth(level, rhs, solution, post_sweeps);
}

} // namespace fem::solver::detail
