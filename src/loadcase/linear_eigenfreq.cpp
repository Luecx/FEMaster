/******************************************************************************
 * @file LinearEigenfrequency.cpp
 * @brief Linear eigenfrequency analysis via affine null-space constraints.
 *
 * Pipeline:
 *   0) Assign sections/materials.
 *   1) Build unconstrained DOF indexing (node×6 → active dof id or -1).
 *   2) Build constraint equations (supports/ties/couplings → C u = d).
 *   3) Assemble active stiffness K (n×n) and mass M (n×n).
 *   4) Build ConstraintTransformer → T, u_p (affine map; u_p=0 for homogeneous).
 *   5) Reduce A = Tᵀ K T,  Mr = Tᵀ M T.
 *   6) Solve A φ = λ Mr φ (smallest λ via shift–invert at σ=0).
 *   7) Recover full modes u = T φ; expand to node×6; compute participations.
 *   8) Print summary & write results.
 *
 * @date    15.09.2025
 * @author  Finn
 ******************************************************************************/

#include "linear_eigenfreq.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

#include "../core/logging.h"
#include "../solve/eigen.h"

#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"

#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"

using fem::constraint::ConstraintTransformer;

namespace fem { namespace loadcase {

//------------------------------------------------------------------------------
// Result container for one mode
//------------------------------------------------------------------------------
struct EigenMode {
    Precision     lambda;        // eigenvalue (ω^2)
    Precision     freq;          // f = sqrt(λ)/(2π)
    DynamicVector q_mode;        // reduced coordinates (q-space)
    DynamicMatrix mode_mat;      // expanded node×6 for writer
    Vec6          participation; // simple modal participation (x,y,z,rx,ry,rz)

    explicit EigenMode(Precision lam, DynamicVector q)
        : lambda(lam), freq(std::sqrt(std::max<Precision>(0, lam)) / (2 * M_PI)), q_mode(std::move(q)) {}
};

//------------------------------------------------------------------------------
// Build 6 axis “unit” direction vectors over active DOFs (length n)
//   col 0..2 → translations x/y/z, col 3..5 → rotations rx/ry/rz
//------------------------------------------------------------------------------
static DynamicMatrix build_active_axis_vectors(const IndexMatrix& active_dof_idx_mat) {
    const int n = active_dof_idx_mat.maxCoeff() + 1; // active system size
    DynamicMatrix U = DynamicMatrix::Zero(n, 6);

    const int nNodes = static_cast<int>(active_dof_idx_mat.rows());
    for (int node = 0; node < nNodes; ++node) {
        for (int dof = 0; dof < 6; ++dof) {
            int gid = active_dof_idx_mat(node, dof);
            if (gid >= 0) {
                // Mark the corresponding axis column
                U(gid, dof) += 1.0;
            }
        }
    }
    return U;
}

// Return an exponent e so that participation values * 10^e are O(1)
static int compute_scaling_exponent(const std::vector<EigenMode>& modes) {
    Precision max_abs = 0;
    for (const auto& m : modes) {
        for (int i = 0; i < 6; ++i)
            max_abs = std::max(max_abs, std::abs(m.participation(i)));
    }
    if (max_abs <= Precision(0)) return 0;
    // choose e = ceil(-log10(max_abs)) so that max_abs * 10^e <= 1
    const double e = std::ceil(-std::log10(static_cast<double>(max_abs)));
    // clamp a bit to avoid silly exponents if values are already order 1
    return static_cast<int>(e);
}

// Pretty, aligned table: eigenvalue, frequency, and (scaled) participations
static void display_eigen_summary(const std::vector<EigenMode>& modes) {
    const int      exp10 = compute_scaling_exponent(modes);
    const Precision scale = std::pow(Precision(10), exp10);

    logging::info(true, "");
    logging::info(true, "Eigenfrequency summary");
    logging::info(true, std::setw(42), "", std::setw(33), "PARTICIPATION", " (x10^", exp10, ")");
    logging::info(true,
                  std::setw(4),  "Idx",
                  std::setw(20), "Eigenvalue",
                  std::setw(18), "Eigenfreq",
                  std::setw(12), "x",  std::setw(8), "y",  std::setw(8), "z",
                  std::setw(8),  "rx", std::setw(8), "ry", std::setw(8), "rz");

    for (size_t i = 0; i < modes.size(); ++i) {
        const auto& m = modes[i];
        logging::info(true,
            std::setw(4),  i + 1,
            std::setw(20), std::scientific, std::setprecision(6), Precision(m.lambda),
            std::setw(18), std::fixed,      std::setprecision(6), Precision(m.freq),
            std::setw(12), std::fixed,      std::setprecision(3), Precision(m.participation(0) * scale),
            std::setw(8),  std::fixed,      std::setprecision(3), Precision(m.participation(1) * scale),
            std::setw(8),  std::fixed,      std::setprecision(3), Precision(m.participation(2) * scale),
            std::setw(8),  std::fixed,      std::setprecision(3), Precision(m.participation(3) * scale),
            std::setw(8),  std::fixed,      std::setprecision(3), Precision(m.participation(4) * scale),
            std::setw(8),  std::fixed,      std::setprecision(3), Precision(m.participation(5) * scale)
        );
    }
    logging::info(true, "");
}

//------------------------------------------------------------------------------
// Write results
//------------------------------------------------------------------------------
static void write_results(const std::vector<EigenMode>& modes,
                          reader::Writer*               writer,
                          int                           loadcase_id)
{
    writer->add_loadcase(loadcase_id);

    DynamicVector eigenvalues(modes.size());
    DynamicVector eigenfreqs (modes.size());

    for (size_t i = 0; i < modes.size(); ++i) {
        eigenvalues(i) = modes[i].lambda;
        eigenfreqs (i) = modes[i].freq;
        writer->write_eigen_matrix(modes[i].mode_mat,      "MODE_SHAPE_"    + std::to_string(i + 1));
        writer->write_eigen_matrix(modes[i].participation, "PARTICIPATION_" + std::to_string(i + 1));
    }
    writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    writer->write_eigen_matrix(DynamicMatrix(eigenfreqs ), "EIGENFREQUENCIES");
}

//------------------------------------------------------------------------------
// LinearEigenfrequency
//------------------------------------------------------------------------------
LinearEigenfrequency::LinearEigenfrequency(ID id,
                                           reader::Writer* writer,
                                           model::Model*   model,
                                           int             numEigenvalues)
    : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

/**
 * @brief Run the linear eigenfrequency analysis (null-space constrained).
 */
void LinearEigenfrequency::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "================================================================================");
    logging::info(true, "LINEAR EIGENFREQUENCY ANALYSIS");
    logging::info(true, "================================================================================");
    logging::info(true, "");

    // (0) Sections/materials
    m_model->assign_sections();

    // (1) Unconstrained system DOF indices (node×6 → active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraints from supports/ties/couplings
    auto equations = Timer::measure(
        [&]() { return m_model->build_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // (3) Assemble active stiffness K and mass M (n×n)
    auto K = Timer::measure(
        [&]() { return m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );
    auto M = Timer::measure(
        [&]() { return m_model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M"
    );

    // (4) Build constraint transformer
    ConstraintTransformer::BuildOptions copt;
    copt.set.scale_columns = true;
    ConstraintTransformer CT(
        equations,
        active_dof_idx_mat, // system DOF map
        K.rows(),           // n
        copt
    );

    // (5) Reduced operators A = Tᵀ K T, Mr = Tᵀ M T
    auto A  = Timer::measure([&]() { return CT.assemble_A(K);  }, "assembling A = T^T K T");
    auto Mr = Timer::measure([&]() { return CT.assemble_B(M);  }, "assembling Mr = T^T M T");

    // (6) Solve generalized EVP A φ = λ Mr φ ; get the *lowest* λ first.
    // Shift-invert at σ=0 makes the smallest λ dominant in SI space.
    solver::EigenOpts eigopt;
    eigopt.mode  = solver::EigenMode::ShiftInvert;
    eigopt.sigma = 0.0;
    eigopt.sort  = solver::EigenOpts::Sort::LargestMagn; // largest in SI ↔ smallest original
    eigopt.ncv   = std::min<int>(int(A.rows()), std::max<int>(3 * num_eigenvalues + 20, num_eigenvalues + 2));

    const int k_req = std::max(1, std::min(num_eigenvalues, int(A.rows())));
    auto eig_pairs = Timer::measure(
        [&]() { return solver::eigs(device, A, Mr, k_req, eigopt); },
        "solving generalized EVP (modal)"
    );

    // Collect modes (ascending λ)
    std::vector<EigenMode> modes; modes.reserve(eig_pairs.size());
    for (const auto& p : eig_pairs) modes.emplace_back(Precision(p.value), p.vector);
    std::sort(modes.begin(), modes.end(),
              [](const EigenMode& a, const EigenMode& b){ return a.lambda < b.lambda; });

    // (7) Expand mode shapes and compute simple participations
    //     u_mode = T q_mode (u_p=0 in modal case)
    DynamicMatrix axes_full = build_active_axis_vectors(active_dof_idx_mat); // n×6
    for (auto& m : modes) {
        // Full vector u
        DynamicVector u_mode = CT.map().apply_T(m.q_mode);
        // Node×6
        m.mode_mat = mattools::expand_vec_to_mat(active_dof_idx_mat, u_mode);
        // Participations p_i = e_iᵀ (M u) using axis columns as e_i
        DynamicVector Mu = M * u_mode;
        for (int i = 0; i < 6; ++i) m.participation(i) = axes_full.col(i).dot(Mu);
    }

    // (8) Print + write
    display_eigen_summary(modes);
    write_results(modes, m_writer, m_id);

    logging::info(true, "Eigenfrequency analysis completed.");
}

}} // namespace fem::loadcase