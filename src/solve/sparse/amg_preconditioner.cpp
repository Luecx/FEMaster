#include "amg_preconditioner.h"

#include "../../core/logging.h"

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

    m_scaling.resize(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        const Precision diagonal = matrix.coeff(i, i);
        logging::error(diagonal > diagonal_epsilon,
                       "AMG requires a positive diagonal, invalid entry at row ", i, ": ", diagonal);
        m_scaling[i] = Precision(1) / std::sqrt(diagonal);
    }

    SparseMatrix scaled = matrix;
    for (Eigen::Index column = 0; column < scaled.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(scaled, column); entry; ++entry) {
            entry.valueRef() *= m_scaling[entry.row()] * m_scaling[entry.col()];
        }
    }
    scaled.makeCompressed();

    m_levels.clear();
    m_workspace.clear();

    Level first;
    first.a = std::move(scaled);
    m_levels.push_back(std::move(first));

    while (static_cast<int>(m_levels.size()) < max_levels
           && m_levels.back().a.rows() > coarse_size) {
        Level& fine = m_levels.back();
        prepare_level(fine);

        int aggregate_count = 0;
        const auto aggregates = build_aggregates(fine.a, fine.inv_diag, aggregate_count);
        if (aggregate_count <= 0 || aggregate_count >= fine.a.rows()) {
            break;
        }

        const SparseMatrix tentative = build_tentative_prolongator(aggregates, aggregate_count);
        const Precision spectral_radius = estimate_spectral_radius(fine.a, fine.inv_diag);
        fine.p = smooth_prolongator(fine.a, tentative, fine.inv_diag, spectral_radius);
        fine.pt = fine.p.transpose();
        fine.pt.makeCompressed();

        SparseMatrix coarse = fine.pt * fine.a * fine.p;
        coarse.makeCompressed();
        if (coarse.rows() >= fine.a.rows()) {
            break;
        }

        Level next;
        next.a = std::move(coarse);
        m_levels.push_back(std::move(next));
    }

    prepare_level(m_levels.back());

    m_workspace.resize(m_levels.size());
    for (std::size_t level = 0; level + 1 < m_levels.size(); ++level) {
        m_workspace[level].residual.resize(m_levels[level].a.rows());
        m_workspace[level].coarse_rhs.resize(m_levels[level + 1].a.rows());
        m_workspace[level].coarse_error.resize(m_levels[level + 1].a.rows());
    }

    m_scaled_rhs.resize(matrix.rows());
    m_scaled_solution.resize(matrix.rows());

    m_coarse_solver.compute(m_levels.back().a);
    logging::error(m_coarse_solver.info() == Eigen::Success,
                   "AMG coarse-level LDLT factorization failed");
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
    DynamicVector inv_diag(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        const Precision diagonal = matrix.coeff(i, i);
        logging::error(diagonal > diagonal_epsilon,
                       "AMG coarse level has a non-positive diagonal at row ", i, ": ", diagonal);
        inv_diag[i] = Precision(1) / diagonal;
    }
    return inv_diag;
}

Precision CpuAmgPreconditioner::estimate_spectral_radius(const SparseMatrix& matrix,
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
    x.normalize();

    Precision rho = 1;
    for (int iteration = 0; iteration < 12; ++iteration) {
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

        int strongest_neighbor = -1;
        Precision strongest_coupling = -std::numeric_limits<Precision>::infinity();

        for (SparseMatrix::InnerIterator entry(matrix, i); entry; ++entry) {
            const int j = static_cast<int>(entry.row());
            if (j == i || aggregates[static_cast<std::size_t>(j)] >= 0) {
                continue;
            }

            const Precision strength = coupling_strength(entry.value(), inv_diag[i], inv_diag[j]);
            if (strength > 0 && strength > strongest_coupling) {
                strongest_coupling = strength;
                strongest_neighbor = j;
            }
        }

        const int aggregate = aggregate_count++;
        aggregates[static_cast<std::size_t>(i)] = aggregate;
        if (strongest_neighbor >= 0) {
            aggregates[static_cast<std::size_t>(strongest_neighbor)] = aggregate;
        }
    }

    return aggregates;
}

SparseMatrix CpuAmgPreconditioner::build_tentative_prolongator(
        const std::vector<int>& aggregates,
        int aggregate_count) {
    TripletList triplets;
    triplets.reserve(aggregates.size());

    for (int row = 0; row < static_cast<int>(aggregates.size()); ++row) {
        triplets.emplace_back(row,
                              aggregates[static_cast<std::size_t>(row)],
                              Precision(1));
    }

    SparseMatrix tentative(static_cast<Eigen::Index>(aggregates.size()), aggregate_count);
    tentative.setFromTriplets(triplets.begin(), triplets.end());
    tentative.makeCompressed();
    return tentative;
}

SparseMatrix CpuAmgPreconditioner::smooth_prolongator(const SparseMatrix& matrix,
                                                       const SparseMatrix& tentative,
                                                       const DynamicVector& inv_diag,
                                                       Precision spectral_radius) {
    SparseMatrix correction = matrix * tentative;
    const Precision omega = Precision(4) / (Precision(3) * spectral_radius);

    for (Eigen::Index column = 0; column < correction.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(correction, column); entry; ++entry) {
            entry.valueRef() *= omega * inv_diag[entry.row()];
        }
    }

    SparseMatrix prolongator = tentative - correction;
    prolongator.makeCompressed();
    return prolongator;
}

void CpuAmgPreconditioner::prepare_level(Level& level) {
    level.inv_diag = inverse_diagonal(level.a);
    const Precision spectral_radius = estimate_spectral_radius(level.a, level.inv_diag);
    level.smoother_omega = Precision(1) / spectral_radius;
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
