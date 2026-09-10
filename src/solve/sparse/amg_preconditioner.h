/**
 * @file amg_preconditioner.h
 * @brief Declares the experimental CPU algebraic multigrid preconditioner.
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
 * @brief Experimental symmetric matrix-only spectral aggregation AMG preconditioner.
 *
 * The hierarchy is constructed once for a fixed SPD matrix. Aggregates are
 * formed only from the scaled matrix graph. On every aggregate a small dense
 * eigenproblem is solved and the locally lowest-energy modes form the coarse
 * basis. No node, element, DOF, or constraint metadata is used.
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
        SparseMatrix p;
        SparseMatrix pt;
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
    static constexpr int pre_sweeps = 1;
    static constexpr int post_sweeps = 1;
    static constexpr int spectral_iterations = 8;

    static constexpr int target_aggregate_size = 12;
    static constexpr int max_aggregation_passes = 6;
    static constexpr int max_local_basis = 6;
    static constexpr Precision local_spectral_cutoff = Precision(0.20);

    static DynamicVector inverse_diagonal(const SparseMatrix& matrix);
    static DynamicVector equilibrate(SparseMatrix& matrix);
    static Precision estimate_spectral_radius(const SparseMatrix& matrix,
                                              const DynamicVector& inv_diag);

    static std::vector<int> build_aggregates(const SparseMatrix& matrix,
                                             int& aggregate_count,
                                             int& min_aggregate_size,
                                             int& max_aggregate_size,
                                             Precision& average_aggregate_size);
    static std::vector<int> pair_graph(const SparseMatrix& graph,
                                       const std::vector<int>& group_sizes,
                                       int& coarse_count,
                                       std::vector<int>& coarse_sizes);
    static SparseMatrix aggregate_strength_graph(const SparseMatrix& graph,
                                                 const std::vector<int>& aggregates,
                                                 int aggregate_count);
    static SparseMatrix build_spectral_prolongator(const SparseMatrix& matrix,
                                                   const std::vector<int>& aggregates,
                                                   int aggregate_count,
                                                   Eigen::Index& coarse_size,
                                                   std::vector<Eigen::Index>& rank_histogram);
    static void scale_prolongator_columns(SparseMatrix& prolongator,
                                          const DynamicVector& scaling);

    void prepare_level(Level& level);
    void smooth(const Level& level,
                const DynamicVector& rhs,
                DynamicVector& solution,
                int sweeps) const;
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
