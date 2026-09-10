#include "amg_preconditioner.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace fem::solver::detail {

namespace {
constexpr Precision diagonal_epsilon = 1e-30;
constexpr int max_aggregate_size = 4;
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

    logging::info(true, "Building CPU aggregation AMG hierarchy");
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
        fine.aggregates = build_aggregates(fine.a, fine.inv_diag, aggregate_count);
        fine.coarse_size = aggregate_count;
        timer.stop();

        logging::info(true,
                      "level ", level_index,
                      " aggregation: ", timer.elapsed(), " ms",
                      " -> N=", aggregate_count);

        if (aggregate_count <= 0 || aggregate_count >= fine.a.rows()) {
            logging::warning(false, "AMG aggregation did not reduce the system; stopping hierarchy construction");
            break;
        }

        timer.start();
        SparseMatrix coarse = build_coarse_matrix(fine.a, fine.aggregates, aggregate_count);
        fine.coarse_scaling = equilibrate(coarse);
        coarse.makeCompressed();
        timer.stop();

        logging::info(true,
                      "level ", level_index,
                      " Galerkin assembly + equilibration: ", timer.elapsed(), " ms",
                      " -> nnz=", coarse.nonZeros());

        if (coarse.rows() >= fine.a.rows()) {
            logging::warning(false, "AMG coarse matrix did not reduce the system; stopping hierarchy construction");
            break;
        }

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

std::vector<int> CpuAmgPreconditioner::build_aggregates(const SparseMatrix& matrix,
                                                         const DynamicVector& inv_diag,
                                                         int& aggregate_count) {
    const int n = static_cast<int>(matrix.rows());
    std::vector<int> aggregates(static_cast<std::size_t>(n), -1);
    aggregate_count = 0;

    for (int i = 0; i < n; ++i) {
        if (aggregates[static_cast<std::size_t>(i)] >= 0) {
            continue;
        }

        std::array<int, max_aggregate_size - 1> neighbors {};
        std::array<Precision, max_aggregate_size - 1> strengths {};
        neighbors.fill(-1);
        strengths.fill(-std::numeric_limits<Precision>::infinity());

        for (SparseMatrix::InnerIterator entry(matrix, i); entry; ++entry) {
            const int j = static_cast<int>(entry.row());
            if (j == i || aggregates[static_cast<std::size_t>(j)] >= 0) {
                continue;
            }

            const Precision strength = coupling_strength(entry.value(), inv_diag[i], inv_diag[j]);
            if (!(strength > 0)) {
                continue;
            }

            for (int slot = 0; slot < max_aggregate_size - 1; ++slot) {
                if (strength <= strengths[static_cast<std::size_t>(slot)]) {
                    continue;
                }
                for (int move = max_aggregate_size - 2; move > slot; --move) {
                    strengths[static_cast<std::size_t>(move)] = strengths[static_cast<std::size_t>(move - 1)];
                    neighbors[static_cast<std::size_t>(move)] = neighbors[static_cast<std::size_t>(move - 1)];
                }
                strengths[static_cast<std::size_t>(slot)] = strength;
                neighbors[static_cast<std::size_t>(slot)] = j;
                break;
            }
        }

        const int aggregate = aggregate_count++;
        aggregates[static_cast<std::size_t>(i)] = aggregate;
        for (const int neighbor : neighbors) {
            if (neighbor >= 0 && aggregates[static_cast<std::size_t>(neighbor)] < 0) {
                aggregates[static_cast<std::size_t>(neighbor)] = aggregate;
            }
        }
    }

    return aggregates;
}

SparseMatrix CpuAmgPreconditioner::build_coarse_matrix(
        const SparseMatrix& matrix,
        const std::vector<int>& aggregates,
        int aggregate_count) {
    const int n = static_cast<int>(matrix.rows());

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
        matrix.nonZeros(), static_cast<Eigen::Index>(aggregate_count) * 96);
    triplets.reserve(static_cast<std::size_t>(reserve_count));

    for (int coarse_column = 0; coarse_column < aggregate_count; ++coarse_column) {
        touched.clear();

        for (int fine_column = first[static_cast<std::size_t>(coarse_column)];
             fine_column >= 0;
             fine_column = next[static_cast<std::size_t>(fine_column)]) {
            for (SparseMatrix::InnerIterator entry(matrix, fine_column); entry; ++entry) {
                const int coarse_row = aggregates[static_cast<std::size_t>(entry.row())];
                if (marker[static_cast<std::size_t>(coarse_row)] != coarse_column) {
                    marker[static_cast<std::size_t>(coarse_row)] = coarse_column;
                    values[static_cast<std::size_t>(coarse_row)] = entry.value();
                    touched.push_back(coarse_row);
                } else {
                    values[static_cast<std::size_t>(coarse_row)] += entry.value();
                }
            }
        }

        std::sort(touched.begin(), touched.end());
        for (const int coarse_row : touched) {
            const Precision value = values[static_cast<std::size_t>(coarse_row)];
            if (value != Precision(0)) {
                triplets.emplace_back(coarse_row, coarse_column, value);
            }
        }
    }

    SparseMatrix coarse(aggregate_count, aggregate_count);
    coarse.setFromTriplets(triplets.begin(), triplets.end());
    coarse.makeCompressed();
    return coarse;
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

void CpuAmgPreconditioner::restrict_residual(const Level& level,
                                             const DynamicVector& fine,
                                             DynamicVector& coarse) const {
    coarse.setZero(level.coarse_size);
    for (Eigen::Index i = 0; i < fine.size(); ++i) {
        const int aggregate = level.aggregates[static_cast<std::size_t>(i)];
        coarse[aggregate] += fine[i];
    }
    coarse.array() *= level.coarse_scaling.array();
}

void CpuAmgPreconditioner::prolongate_and_add(const Level& level,
                                              const DynamicVector& coarse,
                                              DynamicVector& fine) const {
    for (Eigen::Index i = 0; i < fine.size(); ++i) {
        const int aggregate = level.aggregates[static_cast<std::size_t>(i)];
        fine[i] += level.coarse_scaling[aggregate] * coarse[aggregate];
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
    restrict_residual(level, workspace.residual, workspace.coarse_rhs);
    workspace.coarse_error.setZero();

    vcycle(level_index + 1, workspace.coarse_rhs, workspace.coarse_error);

    prolongate_and_add(level, workspace.coarse_error, solution);
    smooth(level, rhs, solution, post_sweeps);
}

} // namespace fem::solver::detail
