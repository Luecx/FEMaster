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
#include "../core/timer.h"
#include "../mattools/append_sparse_block.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"
#include "../solve/eigval/solve_eigval.h"

#include <Eigen/LU>
#include <Eigen/SparseCholesky>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

constexpr Precision pi = Precision(3.141592653589793238462643383279502884L);

template<class Function>
decltype(auto) run_quiet(Function&& function) {
    const bool logging_was_enabled = logging::is_enabled();
    logging::disable();

    try {
        if constexpr (std::is_void_v<std::invoke_result_t<Function>>) {
            std::forward<Function>(function)();
            if (logging_was_enabled) logging::enable();
            return;
        } else {
            auto result = std::forward<Function>(function)();
            if (logging_was_enabled) logging::enable();
            return result;
        }
    } catch (...) {
        if (logging_was_enabled) logging::enable();
        throw;
    }
}

/**
 * Computes the static load-pattern responses only when Kr is positive definite.
 *
 * Modal residual enrichment requires a unique static solution. The generic
 * direct solver falls back to SparseQR for singular matrices and can return
 * finite least-squares displacements that omit rigid-body motion. An explicit
 * Cholesky factorization avoids accepting those invalid Ritz directions.
 *
 * Reuse the same sparse factorization for every load-pattern right-hand side.
 * Returning nullopt means that harmonic response must use the direct solver.
 *
 * @param stiffness Constrained structural stiffness Kr.
 * @param loads Constrained load patterns, one right-hand side per column.
 * @return Static responses or nullopt if the static system is not well posed.
 */
std::optional<DynamicMatrix> solve_static_responses(const SparseMatrix& stiffness,
                                                    const DynamicMatrix& loads) {
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
 * 2. Solve Kr s_j = fr_j for the independently defined load patterns.
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
        // Select natural modes below twice the highest excitation frequency.
        // The generalized eigenvalues are squared angular frequencies.
        const Precision cutoff_frequency  = 2 * *std::max_element(frequencies.begin(), frequencies.end());
        const Precision cutoff_omega      = 2 * pi * cutoff_frequency;
        const Precision cutoff_eigenvalue = cutoff_omega * cutoff_omega;

        logging::error(std::isfinite(cutoff_eigenvalue),
            "LinearHarmonic: modal basis upper eigenvalue must be finite");

        if (cutoff_eigenvalue > 0 && Kr.rows() > 1) {
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
                        device, Kr, Mr, Precision(0), cutoff_eigenvalue, opts, &interval_unresolved);
                },
                "constructing harmonic modal basis");

            logging::info(interval_unresolved,
                "Modal interval exceeds sparse eigenpair coverage; using the direct harmonic solver");

            eigenmode_count = static_cast<Index>(eigenpairs.size());
            reduction_basis.resize(Kr.rows(), eigenmode_count);

            // Normalize eigenvectors with respect to the reduced mass matrix
            for (Index mode = 0; mode < eigenmode_count; ++mode) {
                const DynamicVector& eigenvector = eigenpairs[static_cast<std::size_t>(mode)].vector;
                const Precision      mode_mass   = eigenvector.dot(Mr * eigenvector);

                logging::error(std::isfinite(mode_mass) && mode_mass > 0,
                    "LinearHarmonic: modal mass must be finite and positive");

                reduction_basis.col(mode) = eigenvector / std::sqrt(mode_mass);
            }
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

            // A singular Kr has no unique static responses. Do not allow the
            // general solver's QR fallback to silently discard rigid-body modes.
            const auto static_responses = Timer::measure(
                [&]() { return solve_static_responses(Kr, reduced_loads); },
                "factorizing harmonic static load basis");

            if (!static_responses.has_value()) {
                logging::info(true,
                    "Static stiffness is singular or ill-conditioned; using the direct harmonic solver");
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
            projected_damping   = damping.alpha * projected_mass + damping.beta * projected_stiffness;
        }

        logging::info(true                            , "Modal basis modes          : ", eigenmode_count);
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
        const Precision omega     = 2 * pi * frequency;

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
