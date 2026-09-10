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
 * @brief Experimental symmetric aggregation AMG preconditioner.
 *
 * The hierarchy is constructed once for a fixed SPD matrix. Each application
 * performs one symmetric V-cycle. The implementation intentionally lives on
 * the CPU for now so the hierarchy can be validated independently of CUDA.
 *
 * The transfer operator is piecewise constant. It is stored as a fine-to-coarse
 * aggregate map rather than as sparse P/P^T matrices. This keeps hierarchy
 * construction and V-cycle transfers linear in the matrix/vector sizes and
 * avoids the large temporary fill of generic sparse triple products.
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
        Eigen::Index coarse_size = 0;
        DynamicVector inv_diag;
        Precision smoother_omega = 1;
    };

    struct Workspace {
        mutable DynamicVector residual;
        mutable DynamicVector coarse_rhs;
        mutable DynamicVector coarse_error;
    };

    static constexpr Eigen::Index coarse_direct_size = 512;
    static constexpr int max_levels = 16;
    static constexpr int pre_sweeps = 2;
    static constexpr int post_sweeps = 2;

    static DynamicVector inverse_diagonal(const SparseMatrix& matrix);
    static Precision estimate_jacobi_radius_bound(const SparseMatrix& matrix,
                                                   const DynamicVector& inv_diag);
    static std::vector<int> build_aggregates(const SparseMatrix& matrix,
                                             const DynamicVector& inv_diag,
                                             int& aggregate_count);
    static SparseMatrix build_coarse_matrix(const SparseMatrix& matrix,
                                            const std::vector<int>& aggregates,
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
