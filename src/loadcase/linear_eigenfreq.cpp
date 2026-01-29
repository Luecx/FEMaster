/**
 * @file LinearEigenfrequency.cpp
 * @brief Linear eigenfrequency (modal) analysis via affine null-space constraints.
 *
 * Overview
 * --------
 * We enforce linear constraints C u = d by parameterizing all admissible
 * displacements as u = u_p + T q, where q contains the master DOFs. This keeps
 * the reduced operators symmetric positive definite (for SPD K and feasible
 * constraints) and avoids a saddle-point system with Lagrange multipliers.
 *
 * Pipeline
 * --------
 *   (0) Assign sections/materials in the model.
 *   (1) Build unconstrained DOF indexing (node x 6 -> active dof id or -1).
 *   (2) Build constraint equations from supports/ties/couplings -> C u = d.
 *   (3) Assemble active stiffness K (n x n) and mass M (n x n).
 *   (4) Build ConstraintTransformer -> T and u_p (affine map; u_p = 0 if homogeneous).
 *   (5) Reduce:
 *         A  = T^T K T
 *         Mr = T^T M T           (reduced "mass")
 *   (6) Solve generalized EVP:
 *         A * phi = lambda * Mr * phi
 *       We target the *smallest* eigenvalues by shift-invert at sigma = 0.
 *   (7) Recover full mode shapes:
 *         u_mode = T * phi    (u_p = 0 for homogeneous modal constraints)
 *       Expand to node x 6 and compute simple modal participations.
 *   (8) Print a summary table and write all results.
 *
 * Implementation highlights
 * -------------------------
 * - ConstraintTransformer creation is wrapped in Timer::measure(...) so you
 *   see its cost (which can dominate without zero-column compression).
 * - Reduced operators A and Mr are assembled via ConstraintMap routines. If
 *   you use direct solvers (as you do), explicit A/Mr is appropriate.
 * - The eigensolver is configured for shift-invert at sigma=0 to extract the
 *   lowest frequencies first. Sorting is "largest magnitude" in SI space.
 *
 * Notes
 * -----
 * - The modal case is homogeneous in practice (d = 0), hence u_p = 0, but we
 *   keep the general path and do not assume it; diagnostics will tell you.
 * - The helper "participation" is a simple mass-weighted projection onto the
 *   six global axes over the active DOFs. It is not a full modal effective
 *   mass computation but is often sufficient for quick screening.
 *
 * @date    15.09.2025
 * @author  Finn
 */

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

/**
 * @struct EigenMode
 * @brief Small container that holds one modal eigenpair in reduced space
 *        together with derived/expanded quantities for output.
 *
 * Data members
 * ------------
 * - lambda: eigenvalue of the reduced problem (omega^2).
 * - freq  : physical frequency f = sqrt(lambda) / (2*pi), guarded for lambda < 0.
 * - q_mode: eigenvector in reduced coordinates (size = n_master).
 * - mode_mat: full mode in node x 6 layout, suitable for writer output.
 * - participation: simple mass-weighted projections onto global axes (x,y,z,rx,ry,rz).
 */
struct EigenMode {
    Precision     lambda;        // eigenvalue (omega^2)
    Precision     freq;          // f = sqrt(lambda) / (2 * pi)
    DynamicVector q_mode;        // reduced coordinates (q-space)
    DynamicMatrix mode_mat;      // expanded node x 6 for writer
    Vec6          participation; // simple modal participation (x,y,z,rx,ry,rz)

    explicit EigenMode(Precision lam, DynamicVector q)
        : lambda(lam),
          freq(std::sqrt(std::max<Precision>(0, lam)) / (2 * M_PI)),
          q_mode(std::move(q)) {}
};

/**
 * @brief Build 6 "unit" direction vectors over active DOFs (length n).
 *
 * The returned matrix U (n x 6) has one column per global axis:
 *   - col 0..2 -> translations x/y/z,
 *   - col 3..5 -> rotations    rx/ry/rz.
 *
 * For each active (node, dof) pair, we place a 1.0 in the corresponding axis
 * column at the global active DOF index. This provides a quick basis to
 * compute axis participations by dotting with M*u.
 */
static DynamicMatrix build_active_axis_vectors(const IndexMatrix& active_dof_idx_mat) {
    const int n = active_dof_idx_mat.maxCoeff() + 1; // active system size
    DynamicMatrix U = DynamicMatrix::Zero(n, 6);

    const int nNodes = static_cast<int>(active_dof_idx_mat.rows());
    for (int node = 0; node < nNodes; ++node) {
        for (int dof = 0; dof < 6; ++dof) {
            int gid = active_dof_idx_mat(node, dof);
            if (gid >= 0) {
                U(gid, dof) += 1.0;
            }
        }
    }
    return U;
}

static DynamicMatrix field_to_dynamic(const model::Field& field) {
    DynamicMatrix out(field.rows, field.components);
    for (Index i = 0; i < field.rows; ++i) {
        for (Index j = 0; j < field.components; ++j) {
            out(i, j) = field(i, j);
        }
    }
    return out;
}

/**
 * @brief Compute a power-of-ten scaling exponent so that reported participations
 *        are O(1). If all participations are zero, returns 0.
 */
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

/**
 * @brief Pretty-print a summary table (lambda, frequency, participations).
 */
static void display_eigen_summary(const std::vector<EigenMode>& modes) {
    const int       exp10 = compute_scaling_exponent(modes);
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

/**
 * @brief Write eigenvalues, eigenfrequencies, mode shapes, and participations.
 */
static void write_results(const std::vector<EigenMode>& modes,
                          reader::Writer*               writer,
                          int                           loadcase_id)
{
    writer->add_loadcase(loadcase_id);

    DynamicVector eigenvalues(modes.size());
    DynamicVector eigenfreqs (modes.size());
    DynamicVector freqs      (modes.size());

    for (size_t i = 0; i < modes.size(); ++i) {
        eigenvalues(i) = modes[i].lambda;
        eigenfreqs (i) = modes[i].freq * 2 * M_PI;
        freqs      (i) = modes[i].freq;
        writer->write_eigen_matrix(modes[i].mode_mat,      "MODE_SHAPE_"    + std::to_string(i + 1));
        writer->write_eigen_matrix(modes[i].participation, "PARTICIPATION_" + std::to_string(i + 1));
    }
    writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    writer->write_eigen_matrix(DynamicMatrix(eigenfreqs ), "EIGENFREQUENCIES");
    writer->write_eigen_matrix(DynamicMatrix(freqs      ), "FREQUENCIES");
}

/**
 * @class LinearEigenfrequency
 * @brief Modal analysis entry point (constrained via null-space).
 */
LinearEigenfrequency::LinearEigenfrequency(ID id,
                                           reader::Writer* writer,
                                           model::Model*   model,
                                           int             numEigenvalues)
    : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

/**
 * @brief Execute the modal analysis as described in the file header.
 *
 * Logging
 * -------
 * Mirroring LinearStatic, all heavyweight steps are wrapped with Timer::measure
 * for consistent performance diagnostics in your logs.
 */
void LinearEigenfrequency::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR EIGENFREQUENCY ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // (0) Sections/materials
    model->assign_sections();

    // (1) Unconstrained system DOF indices (node x 6 -> active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraints from supports/ties/couplings
    auto equations = Timer::measure(
        [&]() {
            auto groups = this->model->collect_constraints(active_dof_idx_mat, supps);
            report_constraint_groups(groups);
            return groups.flatten();
        },
        "building constraints"
    );

    // (3) Assemble active stiffness K and mass M (n x n)
    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );
    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M"
    );

    // (4) Build constraint transformer (Set -> Builder -> Map)
    //     Wrapped for timing and to keep the logging style uniform with LinearStatic.
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions copt;
            // Recommended: scale columns of C for robust QR (rank detection).
            copt.set.scale_columns = true;
            // Optional tolerances (keep defaults unless your constraints are ill-conditioned):
            // copt.builder.rank_tol_rel = 1e-12;
            // copt.builder.feas_tol_rel = 1e-10;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,   // system DOF map (node,dof) -> global id or -1
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

    // (5) Reduced operators A = T^T K T, Mr = T^T M T
    // Note: These assemble explicit matrices (direct solver friendly).
    auto A  = Timer::measure(
        [&]() { return CT->assemble_A(K);  },
        "assembling A = T^T K T"
    );
    auto Mr = Timer::measure(
        [&]() { return CT->assemble_B(M);  },
        "assembling Mr = T^T M T"
    );

    // (6) Solve generalized EVP A * phi = lambda * Mr * phi; get the lowest modes.
    // We use shift-invert at sigma = 0.0, so "largest magnitude" in SI space
    // corresponds to the smallest original eigenvalues.
    solver::EigenOpts eigopt;
    eigopt.mode  = solver::EigenMode::ShiftInvert;
    eigopt.sigma = 0.0;
    eigopt.sort  = solver::EigenOpts::Sort::LargestMagn; // largest in SI <-> smallest original

    const int k_req = std::max(1, std::min(num_eigenvalues, int(A.rows())));
    auto eig_pairs = Timer::measure(
        [&]() { return solver::eigs(device, A, Mr, k_req, eigopt); },
        "solving generalized EVP (modal)"
    );

    // Collect modes (ascending lambda)
    std::vector<EigenMode> modes; modes.reserve(eig_pairs.size());
    for (const auto& p : eig_pairs) modes.emplace_back(Precision(p.value), p.vector);
    std::sort(modes.begin(), modes.end(),
              [](const EigenMode& a, const EigenMode& b){ return a.lambda < b.lambda; });

    // (7) Expand mode shapes and compute simple participations
    //     For homogeneous constraints, u_p = 0 and u_mode = T * q_mode.
    DynamicMatrix axes_full = build_active_axis_vectors(active_dof_idx_mat); // n x 6
    for (auto& m : modes) {
        // Full vector u
        DynamicVector u_mode = CT->map().apply_T(m.q_mode);
        CT->post_check_static(K, DynamicVector::Zero(K.rows()), u_mode,
                      /*tol_constraint_rel*/1e-10,
                      /*tol_reduced_rel   */std::numeric_limits<Precision>::infinity(),
                      /*tol_full_rel      */std::numeric_limits<Precision>::infinity());
        // Node x 6 layout for writer
        auto mode_field = mattools::expand_vec_to_mat(active_dof_idx_mat, u_mode);
        m.mode_mat = field_to_dynamic(mode_field);
        // Participations p_i = e_i^T (M u) using axis columns as e_i
        DynamicVector Mu = M * u_mode;
        for (int i = 0; i < 6; ++i) m.participation(i) = axes_full.col(i).dot(Mu);
    }

    // (8) Print + write
    display_eigen_summary(modes);
    write_results(modes, writer, id);

    logging::info(true, "Eigenfrequency analysis completed.");
}

}} // namespace fem::loadcase
