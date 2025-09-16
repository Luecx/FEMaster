/******************************************************************************
 * @file linear_buckling.cpp
 * @brief Linear buckling analysis using affine null-space constraints (u = u_p + T q).
 *
 * Pipeline (null-space path, no Lagrange blocks):
 *   0) Assign sections/materials.
 *   1) Build unconstrained DOF indexing (node×6 → active dof id or -1).
 *   2) Build constraint equations (supports/ties/couplings → C u = d).
 *   3) Assemble global load matrix (node×6) and active stiffness K (n×n).
 *   4) Reduce loads to active RHS f (n×1).
 *   5) Build ConstraintTransformer → get T, u_p (affine map).
 *   6) Pre-buckling static solve in reduced space:
 *        A = Tᵀ K T, b = Tᵀ (f − K u_p), solve A q = b, u = u_p + T q.
 *   7) From u, compute integration-point stresses; assemble geometric K_g (n×n).
 *   8) Reduce K and K_g: A = Tᵀ K T, B = Tᵀ K_g T.
 *   9) Solve generalized EVP: A φ = λ (−B) φ for k modes (shifted if desired).
 *  10) Expand each reduced φ (q-space) → full u-mode → node×6 matrix for output.
 *
 * Remarks
 * -------
 * - The EVP is homogeneous; any inhomogeneous constraints were already accounted
 *   for in the preload via u_p (static step).
 * - We estimate a shift σ from the preload displacement using a Rayleigh quotient
 *   in reduced space to guide the eigensolver (optional but helpful).
 *
 * @date    15.09.2025
 * @author  Finn
 ******************************************************************************/

#include "linear_buckling.h"

#include <algorithm>
#include <iomanip>
#include <limits>

#include "../solve/eigen.h"                 // solve(...), eigs(...)
#include "../core/logging.h"

#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"

#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"

namespace fem { namespace loadcase {

using fem::constraint::ConstraintTransformer;

//------------------------------------------------------------------------------
// Small helpers (local to this TU)
//------------------------------------------------------------------------------
struct BucklingMode {
    Precision     lambda;       ///< Buckling factor (eigenvalue)
    DynamicVector q_mode;       ///< Reduced coordinates (in q-space)
    DynamicMatrix mode_mat;     ///< Expanded (node × DOF) layout for writing

    BucklingMode(Precision lam, DynamicVector q) : lambda(lam), q_mode(std::move(q)) {}
};

/** @brief Cheap Rayleigh quotient estimate λ ≈ (qᵀ A q) / (qᵀ (−B) q). */
static inline double estimate_lambda_rayleigh(const SparseMatrix& A,
                                              const SparseMatrix& B,
                                              const DynamicVector& q) {
    if (q.size() == 0) return 0.0;
    DynamicVector Aq = A * q;
    DynamicVector Bq = B * q;          // we use −B later
    const double num = q.dot(Aq);
    const double den = q.dot(-Bq);
    if (std::abs(den) < 1e-20) return 0.0;
    return num / den;
}

/** @brief Pretty print a short table of buckling factors. */
static void print_buckling_summary(const std::vector<BucklingMode>& modes,
                                   double sigma_used,
                                   int k_requested) {
    logging::info(true, "");
    logging::info(true, "Buckling summary");
    logging::up();
    logging::info(true, "requested modes k : ", k_requested);
    logging::info(true, "shift σ           : ", std::scientific, std::setprecision(3), sigma_used);
    logging::down();

    logging::info(true, std::setw(6), "Idx", std::setw(24), "Buckling factor λ");
    for (size_t i = 0; i < modes.size(); ++i) {
        logging::info(true,
                      std::setw(6), i + 1,
                      std::setw(24), std::fixed, std::setprecision(9),
                      Precision(modes[i].lambda));
    }
}

//------------------------------------------------------------------------------
// LinearBuckling
//------------------------------------------------------------------------------
LinearBuckling::LinearBuckling(ID id,
                               reader::Writer* writer,
                               model::Model* model,
                               int numEigenvalues)
    : LoadCase(id, writer, model)
    , num_eigenvalues(numEigenvalues) {}

/**
 * @brief Execute the linear buckling analysis.
 *
 * See file header for the full step list. This version enforces constraints
 * via the null-space transformer, avoiding a saddle-point system.
 */
void LinearBuckling::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "================================================================================");
    logging::info(true, "LINEAR BUCKLING ANALYSIS");
    logging::info(true, "================================================================================");
    logging::info(true, "");

    // (0) Sections/materials
    m_model->assign_sections();

    // (1) Unconstrained DOF index (node×6 → active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Build constraint equations from supports/ties/couplings
    auto equations = Timer::measure(
        [&]() { return m_model->build_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // (3) Global load matrix (node×6) → keep for reporting if you like
    auto global_load_mat = Timer::measure(
        [&]() { return m_model->build_load_matrix(loads); },
        "building global load matrix"
    );

    // (4) Active stiffness K (n×n)
    auto K = Timer::measure(
        [&]() { return m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );

    // (5) Reduce global loads → active RHS f (n×1)
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix → active RHS vector f"
    );

    // (6) Constraint transformer (Set → Builder → Map)
    ConstraintTransformer::BuildOptions copt;
    copt.set.scale_columns = true; // robust QR
    ConstraintTransformer CT(
        equations,
        active_dof_idx_mat,   // system DOF map
        K.rows(),             // n (active DOFs)
        copt
    );

    // (7) Pre-buckling static solve in reduced space
    auto A = Timer::measure(
        [&]() { return CT.assemble_A(K); },
        "assembling A = T^T K T (preload)"
    );
    auto b = Timer::measure(
        [&]() { return CT.assemble_b(K, f); },
        "assembling b = T^T (f - K u_p) (preload)"
    );

    auto q_pre = Timer::measure(
        [&]() { return solve(device, method, A, b); },
        "solving reduced preload A q = b"
    );
    auto u_pre = Timer::measure(
        [&]() { return CT.recover_u(q_pre); },
        "recovering full preload displacement u"
    );

    // (8) Integration-point stresses from preload → K_g assembly
    IPData ip_stress, ip_strain_unused;
    {
        // Many model APIs expect node-wise displacements (node×6) for IP recovery.
        auto U_mat = Timer::measure(
            [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u_pre); },
            "expanding u to node×DOF for IP stress"
        );
        std::tie(ip_stress, ip_strain_unused) = Timer::measure(
            [&]() { return m_model->compute_ip_stress_strain(U_mat); },
            "computing IP stress/strain for Kg"
        );
    }

    auto Kg = Timer::measure(
        [&]() { return m_model->build_geom_stiffness_matrix(active_dof_idx_mat, ip_stress); },
        "assembling geometric stiffness K_g"
    );

    // (9) Reduced buckling operators
    auto B = Timer::measure(
        [&]() { return CT.assemble_B(Kg); },
        "assembling B = T^T K_g T"
    );

    // Optional: estimate a shift σ from preload (Rayleigh in reduced space)
    double sigma = Timer::measure(
        [&]() { return estimate_lambda_rayleigh(A, B, q_pre) / 1e6; },
        "estimating initial shift σ (Rayleigh)"
    );
    if (sigma <= 0) sigma = - sigma;
    if (!std::isfinite(sigma)) sigma = 1e-6;

    // (10) Solve generalized EVP in reduced space: A φ = λ (−B) φ
    SparseMatrix Bneg = (-1.0) * B;

    const int k_req = std::max(1, std::min<int>(num_eigenvalues, std::max<int>(1, int(A.rows()))));
    solver::EigenOpts eigopt;
    eigopt.mode  = solver::EigenMode::Buckling;
    eigopt.sigma = sigma;
    eigopt.ncv   = std::min<int>(int(A.rows()), std::max<int>(3*k_req + 20, k_req + 2));
    eigopt.sort  = solver::EigenOpts::Sort::LargestMagn;

    auto eig_pairs = Timer::measure(
        [&]() { return solver::eigs(device, A, Bneg, k_req, eigopt); },
        "solving reduced generalized EVP for buckling"
    );

    // Package modes and expand to node×DOF for writing
    std::vector<BucklingMode> modes;
    modes.reserve(eig_pairs.size());
    for (const auto& p : eig_pairs) modes.emplace_back(Precision(p.value), p.vector);

    // Sort ascending by buckling factor (smallest first)
    std::sort(modes.begin(), modes.end(),
      [](const BucklingMode& a, const BucklingMode& b) {
          return a.lambda < b.lambda;
      });


    for (auto& m : modes) {
        // Reduced → full vector with constraints
        DynamicVector u_mode = CT.recover_u(m.q_mode);
        // Full → node×DOF
        m.mode_mat = mattools::expand_vec_to_mat(active_dof_idx_mat, u_mode);
    }

    // Summary
    print_buckling_summary(modes, eigopt.sigma, k_req);

    // (11) Write results
    m_writer->add_loadcase(m_id);
    {
        DynamicVector lambdas(modes.size());
        for (size_t i = 0; i < modes.size(); ++i) {
            lambdas(i) = modes[i].lambda;
            m_writer->write_eigen_matrix(modes[i].mode_mat,
                                         "BUCKLING_MODE_" + std::to_string(i + 1));
        }
        m_writer->write_eigen_matrix(DynamicMatrix(lambdas), "BUCKLING_FACTORS");
    }

    // (12) Small post-checks (optional): projected residual of preload
    {
        DynamicVector r_pre  = K * u_pre - f;           // full active residual (reactions + orthogonal)
        DynamicVector red    = CT.map().apply_Tt(r_pre);
        logging::info(true, "");
        logging::info(true, "Preload post-checks");
        logging::up();
        logging::info(true, "||u_pre||2             : ", u_pre.norm());
        logging::info(true, "||C u_pre - d||2       : ", (CT.set().C * u_pre - CT.set().d).norm());
        logging::info(true, "||K u_pre - f||2       : ", r_pre.norm());
        logging::info(true, "||T^T (K u_pre - f)||2 : ", red.norm());
        logging::down();
    }

    logging::info(true, "Buckling analysis completed.");
}

}} // namespace fem::loadcase
