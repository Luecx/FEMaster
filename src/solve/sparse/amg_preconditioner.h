/**
 * @file amg_preconditioner.h
 * @brief Declares the experimental CPU aggregation AMG preconditioner.
 *
 * @author Finn Eggers
 */

#pragma once

#include "../../core/types_eig.h"

#include <Eigen/SparseCholesky>

#include <cstddef>
#include <vector>

namespace fem::solver::detail {

/**
 * @brief Experimental symmetric algebraic aggregation AMG preconditioner.
 *
 * The hierarchy is constructed once for a fixed SPD matrix. Each application
 * performs one symmetric V-cycle. The implementation intentionally lives on
 * the CPU for now so the hierarchy can be validated independently of CUDA.
 *
 * Fine unknowns are paired purely from the scaled matrix graph. The tentative
 * interpolation keeps the low-energy sign implied by the strongest local
 * coupling instead of assuming that every strongly coupled pair has equal
 * values. No node, element, DOF, or constraint metadata is used.
 */
class CpuAmgPreconditioner {
public:
    void compute(const SparseMatrix& matrix);
    void apply(const DynamicVector& rhs, DynamicVector& solution) const;

    [[nodiscard]] std::size_t level_count() const { return m_levels.size(); }
    [[nodiscard]] Eigen::Index level_size(std::size_t level) const { return m_levels[level].a.rows(); }
    [[nodiscard]] Eigen::Index level_nnz(std::size_t level) const { return m_levels[level].a.nonZeros(); }

private:
    struct Level {
        SparseMatrix a;
        std::vector<int> aggregates;
        DynamicVector interpolation_weights;
        Eigen::Index coarse_size = 0;
        DynamicVector coarse_scaling;
        DynamicVector inv_diag;
        Precision smoother_omega = 1;
    };

    struct Workspace {
        mutable DynamicVector residual;
        mutable DynamicVector coarse_rhs;
        mutable DynamicVector coarse_error;
    };

    static constexpr Eigen::Index coarse_direct_size = 512;
    static constexpr int max_levels = 20;
    static constexpr int pre_sweeps = 2;
    static constexpr int post_sweeps = 2;
    static constexpr int spectral_iterations = 8;

    static DynamicVector inverse_diagonal(const SparseMatrix& matrix);
    static DynamicVector equilibrate(SparseMatrix& matrix);
    static Precision estimate_spectral_radius(const SparseMatrix& matrix,
                                              const DynamicVector& inv_diag);
    static std::vector<int> build_aggregates(const SparseMatrix& matrix,
                                             const DynamicVector& inv_diag,
                                             DynamicVector& interpolation_weights,
                                             int& aggregate_count);
    static SparseMatrix build_coarse_matrix(const SparseMatrix& matrix,
                                            const std::vector<int>& aggregates,
                                            const DynamicVector& interpolation_weights,
                                            int aggregate_count);

    void prepare_level(Level& level);
    void smooth(const Level& level,
                const DynamicVector& rhs,
                DynamicVector& solution,
                int sweeps) const;
    void restrict_residual(const Level& level,
                           const DynamicVector& fine,
                           DynamicVector& coarse) const;
    void prolongate_and_add(const Level& level,
                            const DynamicVector& coarse,
                            DynamicVector& fine) const;
    void vcycle(std::size_t level,
                const DynamicVector& rhs,
                DynamicVector& solution) const;

    std::vector<Level> m_levels;
    mutable std::vector<Workspace> m_workspace;
    Eigen::SimplicialLDLT<SparseMatrix> m_coarse_solver;

    DynamicVector m_scaling;
    mutable DynamicVector m_scaled_rhs;
    mutable DynamicVector m_scaled_solution;
};

} // namespace fem::solver::detail
