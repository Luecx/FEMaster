/******************************************************************************
 * @file linear_buckling.cpp
 * @brief Linear buckling analysis using affine null-space constraints (u = u_p + T q).
 *
 * Overview
 * --------
 * We enforce linear constraints C u = d by parameterizing all admissible
 * displacements as u = u_p + T q, where q are the master DOFs. This keeps the
 * reduced operators symmetric (for SPD K and feasible constraints) and avoids
 * a saddle-point KKT system with Lagrange multipliers.
 *
 * Pipeline (null-space path, no Lagrange blocks)
 * ----------------------------------------------
 *   (0) Assign sections/materials.
 *   (1) Build unconstrained DOF indexing (node x 6 -> active dof id or -1).
 *   (2) Build constraint equations (supports/ties/couplings -> C u = d).
 *   (3) Assemble global load matrix (node x 6) and active stiffness K (n x n).
 *   (4) Reduce loads to active RHS f (n x 1).
 *   (5) Build ConstraintTransformer -> T, u_p (affine map).
 *   (6) Pre-buckling static solve in reduced space:
 *         A = T^T K T,  b = T^T (f - K u_p),  solve A q = b,  u = u_p + T q.
 *   (7) From u, compute IP stresses; assemble geometric stiffness K_g (n x n).
 *   (8) Reduce geometric operator:  B = T^T K_g T  (A is already available).
 *   (9) Solve generalized EVP:  A * phi = lambda * (-B) * phi  for k modes
 *       (optionally with a shift).
 *  (10) Expand each reduced phi (q-space) -> full u-mode -> node x 6 matrix.
 *
 * Remarks
 * -------
 * - The eigenproblem is homogeneous; any inhomogeneous constraints have been
 *   reflected in the pre-buckling preload via u_p.
 * - We estimate a shift sigma from the preload displacement using a Rayleigh
 *   quotient in reduced space (optional but often helpful for convergence).
 * - If you use direct solvers for A and the EVP, building explicit A/B is fine.
 *   For very large systems you could switch to operator mode (matrix-free),
 *   but that is not required here.
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

/******************************************************************************
 * @struct BucklingMode
 * @brief Container for one buckling mode in reduced space with expanded output.
 *
 * Data members
 * ------------
 * - lambda  : buckling factor (eigenvalue of the generalized EVP).
 * - q_mode  : eigenvector in reduced coordinates (size = n_master).
 * - mode_mat: full displacement shape in node x DOF layout for writer output.
 ******************************************************************************/
struct BucklingMode {
    Precision     lambda;       // buckling factor (eigenvalue)
    DynamicVector q_mode;       // reduced coordinates (in q-space)
    DynamicMatrix mode_mat;     // expanded (node x DOF) layout for writing

    BucklingMode(Precision lam, DynamicVector q) : lambda(lam), q_mode(std::move(q)) {}
};

/******************************************************************************
 * @brief Cheap Rayleigh quotient estimate: lambda ~= (q^T A q) / (q^T (-B) q).
 *
 * Purpose
 * -------
 * Provides a quick, scale-informed initial shift estimate from the preload
 * solution (expressed in reduced coordinates). This is not essential, but can
 * improve robustness of the eigensolver in buckling problems.
 ******************************************************************************/
static inline double estimate_lambda_rayleigh(const SparseMatrix& A,
                                              const SparseMatrix& B,
                                              const DynamicVector& q) {
    if (q.size() == 0) return 0.0;
    DynamicVector Aq = A * q;
    DynamicVector Bq = B * q;          // we use -B later
    const double num = q.dot(Aq);
    const double den = q.dot(-Bq);
    if (std::abs(den) < 1e-20) return 0.0;
    return num / den;
}

/******************************************************************************
 * @brief Pretty-print a short table of buckling factors and settings.
 ******************************************************************************/
static void print_buckling_summary(const std::vector<BucklingMode>& modes,
                                   double sigma_used,
                                   int k_requested) {
    logging::info(true, "");
    logging::info(true, "Buckling summary");
    logging::up();
    logging::info(true, "requested modes k : ", k_requested);
    logging::info(true, "shift sigma       : ", std::scientific, std::setprecision(3), sigma_used);
    logging::down();

    logging::info(true, std::setw(6), "Idx", std::setw(24), "Buckling factor lambda");
    for (size_t i = 0; i < modes.size(); ++i) {
        logging::info(true,
                      std::setw(6), i + 1,
                      std::setw(24), std::fixed, std::setprecision(9),
                      Precision(modes[i].lambda));
    }
}

/******************************************************************************
 * @class LinearBuckling
 * @brief Buckling analysis entry point (constrained via null-space).
 ******************************************************************************/
LinearBuckling::LinearBuckling(ID id,
                               reader::Writer* writer,
                               model::Model* model,
                               int numEigenvalues)
    : LoadCase(id, writer, model)
    , num_eigenvalues(numEigenvalues) {}

/******************************************************************************
 * @brief Execute the linear buckling analysis.
 *
 * Implementation notes
 * --------------------
 * - All heavyweight steps are wrapped with Timer::measure for consistent
 *   performance diagnostics (same style as LinearStatic and LinearEigenfrequency).
 * - The ConstraintTransformer creation is wrapped for timing, and diagnostics
 *   (rank, homogeneity, feasibility) are printed.
 ******************************************************************************/
void LinearBuckling::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "================================================================================");
    logging::info(true, "LINEAR BUCKLING ANALYSIS");
    logging::info(true, "================================================================================");
    logging::info(true, "");

    // (0) Sections/materials
    m_model->assign_sections();

    // (1) Unconstrained DOF index (node x 6 -> active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Build constraint equations from supports/ties/couplings
    auto equations = Timer::measure(
        [&]() { return m_model->build_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // (3) Global load matrix (node x 6) -> keep for reporting if you like
    auto global_load_mat = Timer::measure(
        [&]() { return m_model->build_load_matrix(loads); },
        "building global load matrix"
    );

    // (4) Active stiffness K (n x n)
    auto K = Timer::measure(
        [&]() { return m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );

    // (5) Reduce global loads -> active RHS f (n x 1)
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f"
    );

    // (6) Build constraint transformer (Set -> Builder -> Map), wrapped for timing
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions copt;
            // Recommended: scale columns of C for robust QR (rank detection, zero-column compression respected).
            copt.set.scale_columns = true;
            // Optional tolerances (adjust if constraints are ill-conditioned):
            // copt.builder.rank_tol_rel = 1e-12;
            // copt.builder.feas_tol_rel = 1e-10;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,   // system DOF map
                K.rows(),             // n (active DOFs)
                copt
            );
        },
        "building constraint transformer"
    );

    // Diagnostics for the constraint transformer
    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)     : ", CT->report().m);
    logging::info(true, "n (cols of C)     : ", CT->report().n);
    logging::info(true, "rank(C)           : ", CT->rank());
    logging::info(true, "masters (n-r)     : ", CT->n_master());
    logging::info(true, "homogeneous       : ", CT->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", CT->feasible() ? "true" : "false");
    if (!CT->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", CT->report().residual_norm);
    }
    logging::down();

    // (7) Pre-buckling static solve in reduced space
    //     We form A and b explicitly (direct-solver friendly), then solve for q_pre.
    auto A = Timer::measure(
        [&]() { return CT->assemble_A(K); },
        "assembling A = T^T K T (preload)"
    );
    auto b = Timer::measure(
        [&]() { return CT->assemble_b(K, f); },
        "assembling b = T^T (f - K u_p) (preload)"
    );

    auto q_pre = Timer::measure(
        [&]() { return solve(device, method, A, b); },
        "solving reduced preload A q = b"
    );
    auto u_pre = Timer::measure(
        [&]() { return CT->recover_u(q_pre); },
        "recovering full preload displacement u"
    );

    // (8) Integration-point stresses from preload -> K_g assembly
    //     Most model APIs expect node-wise displacements (node x 6) for IP recovery.
    IPData ip_stress, ip_strain_unused;
    {
        auto U_mat = Timer::measure(
            [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u_pre); },
            "expanding u to node x DOF for IP stress"
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

    // (9) Reduced buckling operator B = T^T K_g T
    auto B = Timer::measure(
        [&]() { return CT->assemble_B(Kg); },
        "assembling B = T^T K_g T"
    );

    if (this->sigma != 0) {
        logging::info(true, "");
        logging::info(true, "User-provided shift sigma = ", std::scientific, std::setprecision(3), this->sigma);
    } else {
        // Optional: estimate an initial shift sigma from the preload (Rayleigh in reduced space).
        // We scale down by 1e6 to avoid aggressive shifts on very stiff systems; clamp to a small positive value.
        this->sigma = Timer::measure(
            [&]() { return estimate_lambda_rayleigh(A, B, q_pre) / 1e6; },
            "estimating initial shift sigma (Rayleigh)"
        );
        if (this->sigma <= 0) this->sigma = -this->sigma;     // ensure positive
        if (!std::isfinite(this->sigma)) this->sigma = 1e-6;  // fallback
    }


    // (10) Solve generalized EVP in reduced space: A * phi = lambda * (-B) * phi
    // Many libraries expect the "K x = lambda M x" form; we pass (-B) as the "mass" side.
    SparseMatrix Bneg = (-1.0) * B;

    const int k_req = std::max(1, std::min<int>(num_eigenvalues, std::max<int>(1, int(A.rows()))));
    solver::EigenOpts eigopt;
    eigopt.mode  = solver::EigenMode::Buckling;     // your solver's mode for buckling EVP
    eigopt.sigma = sigma;                           // shift to guide convergence
    eigopt.sort  = solver::EigenOpts::Sort::LargestMagn; // consistent with your eigs(...) contract

    auto eig_pairs = Timer::measure(
        [&]() { return solver::eigs(device, A, Bneg, k_req, eigopt); },
        "solving reduced generalized EVP for buckling"
    );

    // Package modes and expand to node x DOF for writing
    std::vector<BucklingMode> modes;
    modes.reserve(eig_pairs.size());
    for (const auto& p : eig_pairs) modes.emplace_back(Precision(p.value), p.vector);

    // Sort ascending by buckling factor (smallest first is conventional)
    std::sort(modes.begin(), modes.end(),
        [](const BucklingMode& a, const BucklingMode& b) { return a.lambda < b.lambda; });

    for (auto& m : modes) {
        // Reduced -> full vector with constraints
        DynamicVector u_mode = CT->recover_u(m.q_mode);
        // Full -> node x DOF
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
    //      Sanity: residual should mostly live in range(C^T) and T^T r should be small.
    {
        DynamicVector r_pre  = K * u_pre - f;  // full active residual (reactions + orthogonal)
        DynamicVector red    = CT->map().apply_Tt(r_pre);
        logging::info(true, "");
        logging::info(true, "Preload post-checks");
        logging::up();
        logging::info(true, "||u_pre||2             : ", u_pre.norm());
        logging::info(true, "||C u_pre - d||2       : ", (CT->set().C * u_pre - CT->set().d).norm());
        logging::info(true, "||K u_pre - f||2       : ", r_pre.norm());
        logging::info(true, "||T^T (K u_pre - f)||2 : ", red.norm());
        logging::down();
    }

    logging::info(true, "Buckling analysis completed.");
}

}} // namespace fem::loadcase
