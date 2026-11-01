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
 *       We target the *smallest* eigenvalues by shift-invert in physical units.
 *   (7) Recover full mode shapes:
 *         u_mode = T * phi    (u_p = 0 for homogeneous modal constraints)
 *       Expand to node x 6 and project all modes onto global rigid-body motions.
 *   (8) Assemble six mass-centered rigid-body directions and project all modes.
 *   (9) Print a summary table and write all results.
 *
 * Implementation highlights
 * -------------------------
 * - ConstraintTransformer creation is wrapped in Timer::measure(...) so you
 *   see its cost (which can dominate without zero-column compression).
 * - Reduced operators A and Mr are assembled via ConstraintMap routines. If
 *   you use direct solvers (as you do), explicit A/Mr is appropriate.
 * - The negative shift is -0.01 times a quadratic Ritz estimate of the lowest
 *   flexible eigenvalue omega^2, after removing admissible rigid-body modes.
 *   A diagonal-sum scale supplies the shift if the Ritz estimate is unavailable.
 *   Sorting is "largest magnitude" in shift-invert space.
 *
 * Notes
 * -----
 * - The modal case is homogeneous in practice (d = 0), hence u_p = 0, but we
 *   keep the general path and do not assume it; diagnostics will tell you.
 * - Modal participation uses the assembled mass matrix and six rigid-body
 *   translation/rotation directions about the assembled-mass center. All modes
 *   are projected together in expanded DOF space. The reported quantities are
 *   mass-weighted projections, not effective modal masses.
 *
 * @date    15.09.2025
 * @author  Finn
 */

#include "linear_eigenfreq.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>

#include "../core/logging.h"
#include "../solve/eigval/solve_eigval.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/equation.h"

#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/eigenvalue_estimates.h"

using fem::constraint::ConstraintTransformer;

namespace fem { namespace loadcase {
constexpr Precision pi = Precision(3.141592653589793238462643383279502884L);

/**
 * @struct ModalMode
 * @brief Small container that holds one modal eigenpair in reduced space
 *        together with derived/expanded quantities for output.
 *
 * Data members
 * ------------
 * - lambda: eigenvalue of the reduced problem (omega^2).
 * - freq  : physical frequency f = sqrt(lambda) / (2*pi), guarded for lambda < 0.
 * - q_mode: eigenvector in reduced coordinates (size = n_master).
 * - mode_mat: full mode in node x 6 layout, suitable for writer output.
 * - participation: mass-weighted projections onto rigid motions (x,y,z,rx,ry,rz).
 */
struct ModalMode {
    Precision     lambda;        // eigenvalue (omega^2)
    Precision     freq;          // f = sqrt(lambda) / (2 * pi)
    DynamicVector q_mode;        // reduced coordinates (q-space)
    model::Field  mode_mat;      // expanded node x 6 for writer
    Vec6          participation; // simple modal participation (x,y,z,rx,ry,rz)

    explicit ModalMode(Precision lam, DynamicVector q)
        : lambda(lam),
          freq(std::sqrt(std::max<Precision>(0, lam)) / (2 * pi)),
          q_mode(std::move(q)) {}
};

/**
 * Constructs the mass-weighted rigid-body directions for modal participation.
 *
 * Each active node contributes three uniform translations and three rotations.
 * The translational displacement induced by rotation about axis e is
 *
 *     u = e x (position - center).
 *
 * Existing rotational DOFs contribute the corresponding unit rotations.
 * This produces the six excitation directions even for solid-only models,
 * which have no explicit rotational DOFs.
 *
 * The mass center is obtained from the translation/rotation coupling of the
 * assembled mass matrix so that point masses and all assembled inertia are
 * represented. When translational mass is zero, the center is undefined and
 * the global origin is retained. In that case valid pure-rotary modal systems
 * can still be projected without dividing by zero.
 *
 * @param dofs Mapping from nodal components to active global DOF indices.
 * @param positions Global nodal positions.
 * @param M Assembled active-system mass matrix.
 * @return Mass-weighted rigid directions M R with six columns.
 */
static DynamicMatrix build_mass_weighted_rigid_body_directions(const SystemDofIds& dofs,
                                                               const model::Field&  positions,
                                                               const SparseMatrix&  M) {
    // Construct rigid-body translations and rotations about the global origin
    DynamicMatrix R = DynamicMatrix::Zero(M.rows(), 6);

    for (Index node = 0; node < positions.rows; ++node) {
        const Vec3 position = positions.row_vec3(node);
        const Precision x   = position(0);
        const Precision y   = position(1);
        const Precision z   = position(2);

        Mat6 directions = Mat6::Identity();
        directions(0, 4) =  z;
        directions(0, 5) = -y;
        directions(1, 3) = -z;
        directions(1, 5) =  x;
        directions(2, 3) =  y;
        directions(2, 4) = -x;

        // Scatter only existing DOFs; solid nodes have no explicit rotations
        for (int dof = 0; dof < 6; ++dof) {
            const int gid = dofs(node, dof);
            if (gid >= 0) R.row(gid) = directions.row(dof);
        }
    }

    // Determine the center of mass from the assembled inertia, including point
    // masses and translational/rotational coupling already represented in M.
    // For rigid translations, R_t^T M R_t = m I; about the global origin,
    // R_t^T M R_rot = -m [center]x, with [center]x the cross-product matrix.
    DynamicMatrix MR = M * R;

    const Mat3 translational_mass = R.leftCols(3).transpose() * MR.leftCols(3);
    const Mat3 coupling_mass      = R.leftCols(3).transpose() * MR.rightCols(3);
    const Precision total_mass    = translational_mass.trace() / Precision(3);

    logging::error(std::isfinite(total_mass) && total_mass >= Precision(0),
        "Eigenfrequency: invalid total translational mass for participation");

    // With rotary inertia but no translational mass, the center of mass is
    // undefined. Keep the global origin as the rotational reference in this case.
    if (total_mass > Precision(0)) {
        const Precision inv_two_mass = Precision(0.5) / total_mass;
        const Vec3 center(
            (coupling_mass(1, 2) - coupling_mass(2, 1)) * inv_two_mass,
            (coupling_mass(2, 0) - coupling_mass(0, 2)) * inv_two_mass,
            (coupling_mass(0, 1) - coupling_mass(1, 0)) * inv_two_mass
        );

        logging::error(center.allFinite(),
            "Eigenfrequency: invalid center of mass for modal participation");

        // Shift the rotational directions to the mass center. Apply the same
        // transformation directly to M R without another sparse multiplication.
        Mat3 center_cross = Mat3::Zero();
        center_cross(0, 1) = -center(2);
        center_cross(0, 2) =  center(1);
        center_cross(1, 0) =  center(2);
        center_cross(1, 2) = -center(0);
        center_cross(2, 0) = -center(1);
        center_cross(2, 1) =  center(0);

        MR.rightCols(3) += MR.leftCols(3) * center_cross;
    }

    return MR;
}

/**
 * @brief Compute a power-of-ten scaling exponent so that reported participations
 *        are O(1). If all participations are zero, returns 0.
 */
static int compute_scaling_exponent(const std::vector<ModalMode>& modes) {
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
static void display_eigen_summary(const std::vector<ModalMode>& modes) {
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
 * @brief Publish modal frame and step fields through the common output resolver.
 *
 * Eigenvectors are supplied only as DISPLACEMENT. MODE_SHAPE is a derived
 * semantic alias handled by OutputRequestHandler, so the eigensolver does not
 * need a separate mode-shape storage convention.
 */
static void write_results(std::vector<ModalMode>& modes,
                          io::writer::OutputRequestHandler& output,
                          io::writer::ResultWriters* writer,
                          int loadcase_id,
                          model::Model* mdl)
{
    using io::writer::OutputField;

    writer->add_loadcase(loadcase_id, io::writer::WriterStepType::Eigenfrequency);

    const Index num_modes = static_cast<Index>(modes.size());

    model::Field eigenvalues   {"EIGENVALUES",      model::FieldDomain::UNKNOWN, num_modes, 1};
    model::Field eigenfreqs    {"EIGENFREQUENCIES", model::FieldDomain::UNKNOWN, num_modes, 1};
    model::Field freqs         {"FREQUENCIES",      model::FieldDomain::UNKNOWN, num_modes, 1};
    model::Field participation {"PARTICIPATION",    model::FieldDomain::UNKNOWN, num_modes, 6};

    for (Index i = 0; i < num_modes; ++i) {
        const auto index = static_cast<std::size_t>(i);

        eigenvalues(i) = modes[index].lambda;
        eigenfreqs (i) = modes[index].freq * 2 * pi;
        freqs      (i) = modes[index].freq;

        for (Index j = 0; j < 6; ++j) {
            participation(i, j) = modes[index].participation(j);
        }

        output.begin_frame(
            modes[index].freq,
            "_" + std::to_string(i + 1)
        );
        output.provide(OutputField::DISPLACEMENT,        modes[index].mode_mat);
        output.write_frame(*writer, mdl->_data.get());
    }

    // Compact modal metadata belongs to the analysis result itself and is
    // always written, independently of node/element field filtering.
    output.provide_step(OutputField::EIGENVALUES,      eigenvalues);
    output.provide_step(OutputField::EIGENFREQUENCIES, eigenfreqs);
    output.provide_step(OutputField::FREQUENCIES,      freqs);
    output.provide_step(OutputField::PARTICIPATION,    participation);
    output.write_step(*writer, mdl->_data.get());
}

/**
 * @brief Execute the modal analysis as described in the file header.
 *
 * Logging
 * -------
 * Mirroring LinearStatic, all heavyweight steps are wrapped with Timer::measure
 * for consistent performance diagnostics in your logs.
 *
 * Constraint construction reads current direct SUPPORT history and explicitly
 * collector-activated conditions from ModelData. The analysis owns eigenvalue and solver
 * settings only; condition definitions and collector identifiers are not retained.
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
    model->step_begin();

    // (1) Unconstrained system DOF indices (node x 6 -> active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_structural_dof_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraints from supports/ties/couplings
    auto equations = Timer::measure(
        [&]() {
            auto groups = this->model->collect_constraints(active_dof_idx_mat, target_conditions);
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

    // (4) Assemble constraints and build the null-space transformation
    //     Wrapped for timing and to keep the logging style uniform with LinearStatic.
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::Options copt;
            // Optional tolerances (keep defaults unless your constraints are ill-conditioned):
            // copt.null_space.rank_tolerance = 1e-12;
            // copt.null_space.feasibility_tolerance = 1e-10;
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
    logging::info(true, "m (rows of C)     : ", CT->report().equations);
    logging::info(true, "n (cols of C)     : ", CT->report().dofs);
    logging::info(true, "rank(C)           : ", CT->rank());
    logging::info(true, "solver unknowns   : ", CT->unknowns());
    logging::info(true, "homogeneous       : ", CT->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", CT->feasible() ? "true" : "false");
    if (!CT->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", CT->report().residual_norm);
    }
    logging::down();

    // (5) Reduced operators A = T^T K T, Mr = T^T M T
    // Note: These assemble explicit matrices (direct solver friendly).
    auto A  = Timer::measure(
        [&]() { return CT->assemble_system_matrix(K); },
        "assembling A = T^T K T"
    );
    auto Mr = Timer::measure(
        [&]() { return CT->reduce_secondary_matrix(M); },
        "assembling Mr = T^T M T"
    );

    // (6) Solve generalized EVP A * phi = lambda * Mr * phi; get the lowest modes.
    // Estimate the lowest flexible omega^2 from admissible quadratic fields.
    // Matrices, interval bounds and returned eigenvalues remain in physical units.
    CT->initialize_dof_locations(active_dof_idx_mat);
    Precision lambda_estimate = Timer::measure(
        [&]() { return mattools::estimate_eigenvalue_ritz(
            A, Mr, *model->_data->positions, *CT); },
        "estimating lowest flexible omega^2 (quadratic Ritz)"
    );

    // Polynomial samples may vanish at every master, and a point oscillator
    // has no coordinate span at all. Neither implies missing stiffness or mass.
    // Use a geometry-independent characteristic omega^2 in these cases.
    const bool use_diagonal_fallback = !std::isfinite(lambda_estimate) || lambda_estimate <= Precision(0);
    if (use_diagonal_fallback) {
        lambda_estimate = Timer::measure(
            [&]() { return mattools::estimate_spectral_scale(A, Mr); },
            "estimating modal omega^2 scale (diagonal fallback)"
        );
    }

    // A zero scale is valid for free masses with no elastic stiffness. Such
    // systems have no positive physical scale; retain the conventional unit
    // reference for the shift only, leaving all physical eigenvalues unchanged.
    const Precision shift_scale = lambda_estimate == Precision(0) ? Precision(1) : lambda_estimate;

    // Prepare the shift before checking the completed spectral configuration.
    // A strictly negative shift avoids factorizing singular free-body stiffness.
    solver::EigvalOpts eigopt;
    eigopt.mode  = solver::EigvalMode::ShiftInvert;
    eigopt.sigma = -Precision(0.01) * shift_scale;
    eigopt.sort  = solver::EigvalOpts::Sort::LargestMagn;

    // Validate the estimate and shift consecutively before solver/output work.
    logging::error(std::isfinite(lambda_estimate) && lambda_estimate >= Precision(0),
        "Eigenfrequency: neither quadratic Ritz nor diagonal fallback provides a nonnegative eigenvalue scale");
    logging::error(std::isfinite(eigopt.sigma) && eigopt.sigma < 0,
        "Eigenfrequency: invalid negative spectral shift");

    logging::info(use_diagonal_fallback,
        "Quadratic Ritz estimate unavailable; using diagonal-sum fallback.");
    logging::info(lambda_estimate == Precision(0),
        "Zero stiffness scale; using unit reference scale for the negative shift.");
    logging::info(true, use_diagonal_fallback ? "Diagonal omega^2 scale = " : "Quadratic Ritz omega^2 estimate = ",
        lambda_estimate, ", shift sigma = ", eigopt.sigma);

    const int k_req = std::max(1, std::min(num_eigenvalues, int(A.rows())));
    auto eig_pairs = Timer::measure(
        [&]() {
            if (use_eigenvalue_range) {
                return solver::eigvals(device, A, Mr, min_eigenvalue, max_eigenvalue, eigopt);
            }
            return solver::eigvals(device, A, Mr, k_req, eigopt);
        },
        "solving generalized EVP (modal)"
    );

    // Collect the physical eigenvalues and mode shapes (ascending lambda).
    std::vector<ModalMode> modes; modes.reserve(eig_pairs.size());
    for (const auto& p : eig_pairs) modes.emplace_back(Precision(p.value), p.vector);
    std::sort(modes.begin(), modes.end(),
              [](const ModalMode& a, const ModalMode& b) { return a.lambda < b.lambda; });

    // (7) Recover all eigenvectors in the full active DOF space
    DynamicMatrix Phi(K.rows(), static_cast<Index>(modes.size()));
    for (Index i = 0; i < static_cast<Index>(modes.size()); ++i) {
        auto& m = modes[static_cast<std::size_t>(i)];

        // Recover the physical mode and preserve the existing nodal writer output
        DynamicVector u_mode = CT->recover_displacement(m.q_mode);
        CT->post_check_static(K, DynamicVector::Zero(K.rows()), m.q_mode,
                      /*tol_constraint_rel*/1e-10,
                      /*tol_reduced_rel   */std::numeric_limits<Precision>::infinity(),
                      /*tol_full_rel      */std::numeric_limits<Precision>::infinity());

        Phi.col(i) = u_mode;
        m.mode_mat = mattools::expand_vec_to_mat(active_dof_idx_mat, u_mode);
    }

    // (8) Construct the mass-weighted rigid-body directions about the mass center
    const DynamicMatrix MR = build_mass_weighted_rigid_body_directions(
        active_dof_idx_mat, *model->_data->positions, M
    );

    // Project all mass-normalized modes onto the six rigid-body directions
    const DynamicMatrix P = Phi.transpose() * MR;
    for (Index i = 0; i < static_cast<Index>(modes.size()); ++i) {
        modes[static_cast<std::size_t>(i)].participation = P.row(i).transpose();
    }

    // (9) Print + write
    display_eigen_summary(modes);
    write_results(modes, output, writer, id, model);

    logging::info(true, "Eigenfrequency analysis completed.");
    model->step_end();
}
}} // namespace fem::loadcase
