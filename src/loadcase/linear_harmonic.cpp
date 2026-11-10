/**
 * @file linear_harmonic.cpp
 * @brief Implements direct and modal-reduced linear harmonic response analyses.
 *
 * The full constrained system is
 *
 *     (Kr - omega^2 Mr + i omega Cr) q = fr,
 *
 * with Rayleigh damping Cr = alpha Mr + beta Kr. Complex displacements are
 * represented by their real and imaginary parts, so both the direct and the
 * reduced solution use the real block system
 *
 *     [ A -B ] [qr] = [fr]
 *     [ B  A ] [qi]   [ 0],
 *
 * where A = K - omega^2 M and B = omega C in the respective coordinate space.
 *
 * With *MODALBASIS, the solver constructs a mass-orthonormal Ritz space from
 * eigenmodes below twice the maximum excitation frequency, then enriches it
 * with the independent static responses of the harmonic load patterns. All
 * projected matrices remain full: the residual vectors are not eigenmodes.
 * Zero-stiffness modes are identified by their normalized infinity-norm
 * stiffness action. They remain in the dynamic basis, while a mass-orthogonal
 * projection removes them from the static residual problem.
 * The existing full-system direct solver is used if no eigenmode is selected.
 *
 * Load amplitudes and constraints use the step-owned target ConditionState.
 * Direct and collector-activated loads are evaluated at every frequency,
 * whereas the stiffness, mass, damping and constraint mapping stay unchanged.
 *
 * @see LinearHarmonic
 * @see src/loadcase/linear_harmonic.h
 * @author Finn Eggers
 * @date 05.08.2026
 */

#include "linear_harmonic.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../core/logging.h"
#include "../core/run_quiet.h"
#include "../core/timer.h"
#include "../mattools/append_sparse_block.h"
#include "../mattools/projection.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../math/pi.h"
#include "../model/model.h"
#include "../solve/eigval/solve_eigval.h"

#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SparseCholesky>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

/**
 * Computes static load-pattern responses outside identified zero-stiffness modes.
 *
 * Without zero modes, the existing positive-definite sparse solve is retained.
 * With zero modes R, form the mass-orthogonal projector P and solve
 *
 *     (P^T Kr P) z = P^T fr,    s = P z.
 *
 * The projected stiffness still has nullity dim(R). A rank-revealing complete
 * orthogonal decomposition solves for all load patterns without inverting its
 * zero pivots. Applying P to the solution fixes R^T Mr s = 0. An unexpected
 * rank or excessive projected equilibrium residual requests the direct sweep.
 * The physical zero-mode response remains in the dynamic modal basis.
 *
 * @param stiffness Constrained structural stiffness Kr.
 * @param mass Constrained mass matrix Mr in the same coordinates.
 * @param loads Constrained load patterns, one right-hand side per column.
 * @param zero_modes Independent positive-mass modes identified by Kr phi = 0.
 * @return Static responses or nullopt if the complement cannot be resolved.
 */
std::optional<DynamicMatrix> solve_static_responses(const SparseMatrix& stiffness,
                                                    const SparseMatrix& mass,
                                                    const DynamicMatrix& loads,
                                                    const DynamicMatrix& zero_modes) {
    if (zero_modes.cols() > 0) {
        // No static enrichment is needed when there are no load patterns.
        if (loads.cols() == 0) return DynamicMatrix(stiffness.rows(), 0);

        // Assemble the requested projected system explicitly. P acts on
        // displacement coordinates; virtual work projects forces with P^T.
        const DynamicMatrix projector = mattools::build_projection_matrix(zero_modes, mass);
        const DynamicMatrix projected_stiffness = projector.transpose() * (stiffness * projector);
        const DynamicMatrix projected_loads     = projector.transpose() * loads;
        if (!projected_stiffness.allFinite() || !projected_loads.allFinite()) return std::nullopt;

        // The projector preserves the matrix dimension and introduces exactly
        // dim(R) null directions. Extra null directions indicate incomplete
        // modal coverage or an unresolved elastic complement; do not silently
        // accept a least-squares response that loses additional motion.
        Eigen::CompleteOrthogonalDecomposition<DynamicMatrix> factorization;
        factorization.setThreshold(Precision(1e-10));
        factorization.compute(projected_stiffness);
        if (factorization.info() != Eigen::Success
            || factorization.rank() != stiffness.rows() - zero_modes.cols()) return std::nullopt;

        const DynamicMatrix coordinates = factorization.solve(projected_loads);
        const DynamicMatrix responses   = projector * coordinates;
        if (!coordinates.allFinite() || !responses.allFinite()) return std::nullopt;

        // Check equilibrium in the projected test space and the chosen
        // mass-orthogonal displacement gauge. The original loads provide the
        // force scale even when their projected elastic part is nearly zero.
        const DynamicMatrix equilibrium = projector.transpose() * (stiffness * responses - loads);
        const DynamicMatrix weighted    = mass * responses;
        const DynamicMatrix zero_motion = zero_modes.transpose() * weighted;
        const Precision     tolerance   = Precision(1e-8);
        if (!equilibrium.allFinite() || !weighted.allFinite() || !zero_motion.allFinite()) return std::nullopt;
        if (equilibrium.norm() > tolerance * loads.norm()) return std::nullopt;
        if (zero_motion.norm() > tolerance * zero_modes.norm() * weighted.norm()) return std::nullopt;

        return responses;
    }

    // Require a positive-definite constrained stiffness (no rigid-body modes)
#ifdef USE_MKL
    Eigen::PardisoLLT<SparseMatrix> factorization;
#else
    Eigen::SimplicialLLT<SparseMatrix> factorization;
#endif
    factorization.compute(stiffness);
    if (factorization.info() != Eigen::Success) return std::nullopt;

    // The factorization is still needed for zero load patterns to detect
    // rigid-body modes, but there is no right-hand side to solve.
    if (loads.cols() == 0) return DynamicMatrix(stiffness.rows(), 0);

    // Solve all static responses with one factorization, without QR fallback
    const DynamicMatrix responses = factorization.solve(loads);
    if (factorization.info() != Eigen::Success || !responses.allFinite()) {
        return std::nullopt;
    }

    // Reject numerically invalid static solutions before they enter the basis
    const Precision load_norm     = loads.norm();
    const Precision residual_norm = (stiffness * responses - loads).norm();
    const Precision tolerance     = Precision(1e-8);

    if (!std::isfinite(residual_norm)) return std::nullopt;
    if (residual_norm > tolerance * load_norm) return std::nullopt;

    return responses;
}

/**
 * Builds the real representation [A,-B;B,A] of the complex operator A+iB.
 *
 * Both inputs must be square and equally sized. The first block coordinate
 * holds real components and the second holds imaginary components. Scaled
 * sparse blocks are collected in a shared triplet buffer before constructing
 * an independent compressed matrix; neither input is modified.
 *
 * @param A Real part of the complex operator.
 * @param B Imaginary part of the complex operator.
 * @return Real sparse operator with twice the input row and column count.
 */
SparseMatrix build_block_matrix(const SparseMatrix& A, const SparseMatrix& B) {
    logging::error(A.rows() == A.cols(),
        "LinearHarmonic: A must be square");
    logging::error(B.rows() == B.cols(),
        "LinearHarmonic: B must be square");
    logging::error(A.rows() == B.rows(),
        "LinearHarmonic: A and B size mismatch");

    // Reserve all four blocks before embedding real and imaginary coefficients.
    const Index n = static_cast<Index>(A.rows());
    TripletList triplets;
    triplets.reserve(2 * (A.nonZeros() + B.nonZeros()));

    // (A+iB)(u_r+i u_i) gives A*u_r-B*u_i and B*u_r+A*u_i.
    mattools::append_sparse_block(triplets, A, 0, 0);
    mattools::append_sparse_block(triplets, A, n, n);
    mattools::append_sparse_block(triplets, B, 0, n, Precision(-1));
    mattools::append_sparse_block(triplets, B, n, 0);

    SparseMatrix result(2 * n, 2 * n);
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

DynamicVector build_block_rhs(const DynamicVector& real, const DynamicVector& imag) {
    logging::error(real.size() == imag.size(),
        "LinearHarmonic: RHS size mismatch");

    DynamicVector result(2 * real.size());
    result.head(real.size()) = real;
    result.tail(real.size()) = imag;
    return result;
}

} // anonymous namespace

/**
 * Solves the frequency-domain response for every requested excitation frequency.
 *
 * The structural stiffness and lumped mass are assembled once and transformed
 * into homogeneous constraint coordinates. Rayleigh damping is constructed in
 * the same coordinates. If *MODALBASIS is enabled, a single reduced Ritz space
 * is assembled before the sweep:
 *
 * 1. Compute eigenvectors of (Kr, Mr) below the specified cutoff and scale
 *    each vector to unit mass, phi^T Mr phi = 1.
 * 2. Identify zero-stiffness eigenvectors using
 *    ||Kr phi||_inf / (||Kr||_inf ||phi||_inf), retaining them dynamically.
 *    Solve the static load patterns in their mass-orthogonal complement;
 *    without zero modes, retain the ordinary Kr s_j = fr_j solve.
 * 3. Remove the represented components from every s_j using two passes of
 *    Mr-orthogonalization; retain only independent residual directions.
 * 4. Project Kr, Mr and Cr onto the combined basis. In particular, the
 *    projected stiffness and mass need not be diagonal.
 *
 * For each frequency the complex dynamic system is solved either in the
 * reduced Ritz space or directly in constrained coordinates. Real and imaginary
 * displacements are then transformed back to global nodal fields for output.
 * Target conditions and load amplitudes are evaluated at each frequency.
 *
 * @note The Ritz solution is an approximation of the full harmonic response.
 *       If the modal basis is empty, the existing direct solver is used.
 */
void LinearHarmonic::run() {
    // Apply step-specific defaults only when no output variables were requested
    if (output.empty()) {
        using io::writer::OutputField;
        output.request({
            {OutputField::DISPLACEMENT_REAL, model::FieldDomain::NODE},
            {OutputField::DISPLACEMENT_IMAG, model::FieldDomain::NODE},
            {OutputField::STRESS_REAL      , model::FieldDomain::NODE},
            {OutputField::STRESS_IMAG      , model::FieldDomain::NODE},
            {OutputField::STRAIN_REAL      , model::FieldDomain::NODE},
            {OutputField::STRAIN_IMAG      , model::FieldDomain::NODE}
        });
    }

    // Validate the requested analysis procedure and excitation frequencies
    logging::error(!frequencies.empty(),
        "LinearHarmonic: no frequencies defined");
    logging::error(std::all_of(frequencies.begin(), frequencies.end(), [](Precision frequency) {
            return std::isfinite(frequency) && frequency >= 0.0;
        }),
        "LinearHarmonic: frequencies must be finite and non-negative");
    logging::error(constraint_method == ConstraintTransformer::Method::NullSpace,
        "LinearHarmonic: only NULLSPACE constraints are currently supported");
    logging::error(method == solver::DIRECT,
        "LinearHarmonic: only DIRECT solver method is currently supported");

    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR HARMONIC RESPONSE ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Initialize structural sections and the state used by the harmonic step
    model->assign_sections();
    model->step_begin();

    // Enumerate active structural DOFs and collect their support equations
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_structural_dof_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat, target_conditions); },
        "building constraints");

    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // Assemble the real stiffness and mass in the original active DOF space
    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K");

    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M");

    // Construct the homogeneous null-space mapping from reduced to active DOFs
    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations, active_dof_idx_mat, K.rows(), options);
        },
        "building constraint transformer");

    logging::error(transformer->homogeneous(),
        "LinearHarmonic: harmonic constraints must be homogeneous");
    logging::error(transformer->feasible(),
        "LinearHarmonic: constraint system is infeasible");

    // Reduce K and M once, then construct Rayleigh damping Cr = alpha Mr + beta Kr
    auto Kr = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(K); },
        "assembling reduced stiffness matrix Kr");

    auto Mr = Timer::measure(
        [&]() { return transformer->reduce_secondary_matrix(M); },
        "assembling reduced mass matrix Mr");

    auto Cr = Timer::measure(
        [&]() { return damping.build(Mr, Kr); },
        "constructing reduced Rayleigh damping matrix Cr");

    // Keep the Ritz basis separate from the reduced physical matrices:
    // its columns contain both eigenmodes and static residual directions.
    DynamicMatrix reduction_basis(Kr.rows(), 0);
    DynamicMatrix projected_stiffness;
    DynamicMatrix projected_mass;
    DynamicMatrix projected_damping;
    Index         eigenmode_count = 0;

    if (modal_basis) {
        DynamicMatrix      zero_modes(Kr.rows(), 0);
        std::vector<Index> zero_mode_indices;

        // Select natural modes below twice the highest excitation frequency.
        // The generalized eigenvalues are squared angular frequencies.
        const Precision cutoff_frequency  = 2 * *std::max_element(frequencies.begin(), frequencies.end());
        const Precision cutoff_omega      = 2 * math::pi * cutoff_frequency;
        const Precision cutoff_eigenvalue = cutoff_omega * cutoff_omega;

        logging::error(std::isfinite(cutoff_eigenvalue),
            "LinearHarmonic: modal basis upper eigenvalue must be finite");

        if (cutoff_eigenvalue > 0 && Kr.rows() > 1) {
            // The interval is open. Include zero and small negative numerical
            // eigenvalues rather than discarding them before the Ku test.
            const Precision minimum_eigenvalue = std::nextafter(
                -Precision(1e-4) * cutoff_eigenvalue, -std::numeric_limits<Precision>::infinity());
            solver::EigvalOpts opts;
            opts.mode = solver::EigvalMode::ShiftInvert;
            opts.sort = solver::EigvalOpts::Sort::LargestMagn;

            // The sparse eigensolver cannot enumerate all N modes. If its
            // interval includes the complete spectrum (or cannot be certified
            // using at most N-1 eigenpairs), use the exact direct sweep instead.
            bool interval_unresolved = false;
            auto eigenpairs = Timer::measure(
                [&]() {
                    return solver::eigvals(
                        device, Kr, Mr, minimum_eigenvalue, cutoff_eigenvalue, opts, &interval_unresolved);
                },
                "constructing harmonic modal basis");

            logging::info(interval_unresolved,
                "Modal interval exceeds sparse eigenpair coverage; using the direct harmonic solver");

            eigenmode_count = static_cast<Index>(eigenpairs.size());
            reduction_basis.resize(Kr.rows(), eigenmode_count);
            zero_modes.resize(Kr.rows(), eigenmode_count);

            // The induced matrix infinity norm is the maximum absolute row
            // sum, not the largest entry. Assemble row sums explicitly because
            // Kr uses column-major sparse storage. No backend inertia is used
            // to classify or count the returned zero-stiffness directions.
            DynamicVector stiffness_row_sums = DynamicVector::Zero(Kr.rows());
            for (Eigen::Index col = 0; col < Kr.outerSize(); ++col) {
                for (SparseMatrix::InnerIterator entry(Kr, col); entry; ++entry) {
                    stiffness_row_sums(entry.row()) += std::abs(entry.value());
                }
            }
            const Precision stiffness_norm      = stiffness_row_sums.maxCoeff();
            const Precision zero_mode_tolerance = Precision(1e-10);
            logging::error(stiffness_row_sums.allFinite() && std::isfinite(stiffness_norm),
                "LinearHarmonic: stiffness infinity norm must be finite");

            // Normalize eigenvectors with respect to the reduced mass matrix
            for (Index mode = 0; mode < eigenmode_count; ++mode) {
                const DynamicVector& eigenvector = eigenpairs[static_cast<std::size_t>(mode)].vector;
                const Precision      mode_mass   = eigenvector.dot(Mr * eigenvector);

                logging::error(std::isfinite(mode_mass) && mode_mass > 0,
                    "LinearHarmonic: modal mass must be finite and positive");

                reduction_basis.col(mode) = eigenvector / std::sqrt(mode_mass);

                // Test eta = ||Kr phi||_inf / (||Kr||_inf ||phi||_inf).
                // Scale phi first to avoid overflow in the norm product.
                // For Kr = 0 the equivalent inequality also accepts exact
                // zero actions without dividing by a zero matrix norm.
                const Precision mode_norm = reduction_basis.col(mode).lpNorm<Eigen::Infinity>();
                logging::error(std::isfinite(mode_norm) && mode_norm > 0,
                    "LinearHarmonic: invalid mode in zero-stiffness classification");

                const DynamicVector action = Kr * (reduction_basis.col(mode) / mode_norm);
                logging::error(action.allFinite(),
                    "LinearHarmonic: zero-stiffness action must be finite");

                if (action.lpNorm<Eigen::Infinity>() <= zero_mode_tolerance * stiffness_norm) {
                    zero_modes.col(static_cast<Eigen::Index>(zero_mode_indices.size())) = reduction_basis.col(mode);
                    zero_mode_indices.push_back(mode);
                }
            }
            zero_modes.conservativeResize(Eigen::NoChange, static_cast<Eigen::Index>(zero_mode_indices.size()));
        }

        if (eigenmode_count > 0) {
            // Gather independent target load patterns before evaluating their
            // frequency-dependent amplitudes. Each pattern supplies one RHS.
            const auto load_patterns = model->build_load_basis(target_conditions);
            const Index pattern_count = static_cast<Index>(load_patterns.size());
            DynamicMatrix reduced_loads(Kr.rows(), pattern_count);

            for (Index pattern = 0; pattern < pattern_count; ++pattern) {
                const auto& source = load_patterns[static_cast<std::size_t>(pattern)].second;
                const DynamicVector force = mattools::reduce_mat_to_vec(active_dof_idx_mat, source);

                reduced_loads.col(pattern) = transformer->assemble_system_rhs(K, force);
            }

            // Retain zero modes dynamically and project only the static
            // enrichment problem. Unidentified singularities still trigger
            // the existing direct harmonic fallback.
            const auto static_responses = Timer::measure(
                [&]() { return solve_static_responses(Kr, Mr, reduced_loads, zero_modes); },
                "factorizing harmonic static load basis");

            if (!static_responses.has_value()) {
                logging::info(true,
                    "Static stiffness complement is unresolved or ill-conditioned; using the direct harmonic solver");
                reduction_basis.resize(Kr.rows(), 0);
                eigenmode_count = 0;
            } else {
                // Remove the eigenmode and earlier residual contributions using
                // two passes of mass-orthogonalization for numerical stability.
                const Precision residual_tolerance = Precision(1e-10);
                for (Index pattern = 0; pattern < pattern_count; ++pattern) {
                    DynamicVector residual = static_responses->col(pattern);
                    const Precision static_mass = residual.dot(Mr * residual);

                    logging::error(std::isfinite(static_mass) && static_mass >= 0,
                        "LinearHarmonic: static response mass norm must be finite and non-negative");
                    if (static_mass == 0) continue;

                    for (int pass = 0; pass < 2; ++pass) {
                        const DynamicVector coefficients = reduction_basis.transpose() * (Mr * residual);
                        residual -= reduction_basis * coefficients;
                    }

                    // Keep only independent static directions with appreciable
                    // mass norm relative to their unprojected static response.
                    const Precision residual_mass = residual.dot(Mr * residual);

                    logging::error(std::isfinite(residual_mass) && residual_mass >= 0,
                        "LinearHarmonic: residual mass norm must be finite and non-negative");
                    if (residual_mass <= residual_tolerance * residual_tolerance * static_mass) continue;

                    const Eigen::Index column = reduction_basis.cols();
                    reduction_basis.conservativeResize(Kr.rows(), column + 1);
                    reduction_basis.col(column) = residual / std::sqrt(residual_mass);
                }
            }
        }

        if (reduction_basis.cols() > 0) {
            // Static residual vectors are not eigenvectors: project full
            // operators rather than assuming diagonal modal coefficients.
            projected_stiffness = reduction_basis.transpose() * (Kr * reduction_basis);
            projected_mass      = reduction_basis.transpose() * (Mr * reduction_basis);

            // Identified numerical zeros carry inertia and mass damping but
            // no elastic stiffness. Remove their roundoff stiffness entries;
            // at omega = 0 a free zero mode remains physically singular.
            for (Index mode : zero_mode_indices) {
                projected_stiffness.row(mode).setZero();
                projected_stiffness.col(mode).setZero();
            }
            projected_damping   = damping.alpha * projected_mass + damping.beta * projected_stiffness;
        }

        logging::info(true                            , "Modal basis modes          : ", eigenmode_count);
        logging::info(true                            , "Detected zero modes        : ", zero_modes.cols());
        logging::info(true                            , "Residual basis vectors     : ", reduction_basis.cols() - eigenmode_count);
        logging::info(true                            , "Modal basis upper frequency: ", cutoff_frequency);
        logging::info(reduction_basis.cols() == 0     , "Modal basis is empty; using the direct harmonic solver");
    }

    // Write each harmonic response as a dynamic result frame
    writer->add_loadcase(id, io::writer::WriterStepType::Dynamic);

    logging::info(true, "Frequency sweep");
    logging::info(true, std::setw(6), "Idx", std::setw(16), "Frequency", std::setw(18), "Response norm");
    logging::info(true, "----------------------------------------");

    Timer sweep_timer;
    sweep_timer.start();

    // The matrices and Ritz space remain fixed throughout the entire sweep;
    // only the frequency-dependent load and dynamic stiffness change.
    for (Index i = 0; i < static_cast<Index>(frequencies.size()); ++i) {
        const Precision frequency = frequencies[static_cast<std::size_t>(i)];
        const Precision omega     = 2 * math::pi * frequency;

        DynamicVector u_real;
        DynamicVector u_imag;

        run_quiet([&]() {
            // Evaluate direct and collector-activated target loads at the
            // current frequency, then reduce forces into constraint coordinates.
            model::Field global_load_mat = build_target_load_matrix(frequency);
            const DynamicVector force = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat);
            const DynamicVector reduced_force = transformer->assemble_system_rhs(K, force);

            DynamicVector q_real;
            DynamicVector q_imag;

            if (reduction_basis.cols() > 0) {
                // Galerkin projection produces a complex system of dimension r,
                // represented here as a real 2r x 2r dense block matrix.
                const Index basis_size = static_cast<Index>(reduction_basis.cols());
                const DynamicMatrix A  = projected_stiffness - omega * omega * projected_mass;
                const DynamicMatrix B  = omega * projected_damping;

                DynamicMatrix block_matrix(2 * basis_size, 2 * basis_size);
                block_matrix.topLeftCorner    (basis_size, basis_size) =  A;
                block_matrix.topRightCorner   (basis_size, basis_size) = -B;
                block_matrix.bottomLeftCorner (basis_size, basis_size) =  B;
                block_matrix.bottomRightCorner(basis_size, basis_size) =  A;

                // Harmonic forcing is real, so the imaginary right-hand side
                // vanishes. Project the real force onto the Ritz space.
                DynamicVector rhs = DynamicVector::Zero(2 * basis_size);
                rhs.head(basis_size) = reduction_basis.transpose() * reduced_force;

                // Solve the small dense problem and reconstruct both real and
                // imaginary displacements in constrained physical coordinates.
                Eigen::FullPivLU<DynamicMatrix> reduced_solver(block_matrix);
                logging::error(reduced_solver.isInvertible(),
                    "LinearHarmonic: projected dynamic stiffness is singular");

                const DynamicVector solution = reduced_solver.solve(rhs);
                logging::error(solution.allFinite(),
                    "LinearHarmonic: projected harmonic solution contains invalid values");

                q_real = reduction_basis * solution.head(basis_size);
                q_imag = reduction_basis * solution.tail(basis_size);
            } else {
                // No modal basis: solve the full constrained dynamic stiffness
                // using the existing general sparse real block solver.
                SparseMatrix A = Kr - omega * omega * Mr;
                SparseMatrix B = omega * Cr;
                A.makeCompressed();
                B.makeCompressed();

                SparseMatrix block_matrix = build_block_matrix(A, B);
                const DynamicVector zero        = DynamicVector::Zero(reduced_force.size());
                const DynamicVector rhs         = build_block_rhs(reduced_force, zero);

                const DynamicVector solution = solve(
                    device, method, block_matrix, rhs, solver::DirectSolverMatrixType::General);

                const Index reduced_size = static_cast<Index>(reduced_force.size());
                q_real = solution.head(reduced_size);
                q_imag = solution.tail(reduced_size);
            }

            // Recover active DOFs from null-space coordinates and expand the
            // real and imaginary parts into separate global nodal fields.
            u_real = transformer->recover_displacement(q_real);
            u_imag = transformer->recover_displacement(q_imag);

            auto displacement_real = mattools::expand_vec_to_mat(active_dof_idx_mat, u_real);
            auto displacement_imag = mattools::expand_vec_to_mat(active_dof_idx_mat, u_imag);

            // Provide the solved displacements and current loading for lazy
            // recovery of harmonic stress/strain and other requested results.
            using io::writer::OutputField;
            output.begin_frame(frequency, "_" + std::to_string(i + 1));
            output.provide(OutputField::DISPLACEMENT_REAL,   displacement_real);
            output.provide(OutputField::DISPLACEMENT_IMAG,   displacement_imag);
            output.provide(OutputField::EXTERNAL_FORCES,     global_load_mat);
            output.write_frame(*writer, model->_data.get());
        });

        // Report the Euclidean norm of both displacement components together
        logging::info(true, std::setw(6), i + 1,
                      std::setw(16), std::fixed, std::setprecision(6), frequency,
                      std::setw(18), std::scientific, std::setprecision(6),
                      std::sqrt(u_real.squaredNorm() + u_imag.squaredNorm()));
    }

    sweep_timer.stop();
    logging::info(true, "Harmonic frequency sweep elapsed: ", sweep_timer.elapsed(), " ms");

    // Finish the structural step after all frequencies have been written
    logging::info(true, "");
    model->step_end();
}

} // namespace loadcase
} // namespace fem
