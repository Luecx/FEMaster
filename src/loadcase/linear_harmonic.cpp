/**
 * @file linear_harmonic.cpp
 * @brief Implements the direct linear harmonic response load case.
 *
 * The complex system
 *
 *     (A + i B) (u_r + i u_i) = f_r + i f_i
 *
 * is solved as the real block system
 *
 *     [ A -B ] [u_r] = [f_r]
 *     [ B  A ] [u_i]   [f_i]
 *
 * with A = K - omega^2 M and B = omega C. The implementation assumes real
 * reference loads, hence f_i = 0, and Rayleigh damping C = alpha M + beta K.
 * Named load amplitudes are evaluated at the current sweep frequency so the same
 * scalar `Amplitude` representation can describe frequency-dependent forcing.
 *
 * @see src/loadcase/linear_harmonic.h
 * @author Finn Eggers
 * @date 05.08.2026
 */

#include "linear_harmonic.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"
#include "../solve/eigval/solve_eigval.h"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <iomanip>
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

SparseMatrix build_block_matrix(const SparseMatrix& A, const SparseMatrix& B) {
    logging::error(A.rows() == A.cols(),
        "LinearHarmonic: A must be square");
    logging::error(B.rows() == B.cols(),
        "LinearHarmonic: B must be square");
    logging::error(A.rows() == B.rows(),
        "LinearHarmonic: A and B size mismatch");

    using Triplet = Eigen::Triplet<Precision, Index>;

    const Index n = static_cast<Index>(A.rows());
    std::vector<Triplet> triplets;
    triplets.reserve(2 * (A.nonZeros() + B.nonZeros()));

    for (Index outer = 0; outer < static_cast<Index>(A.outerSize()); ++outer) {
        for (SparseMatrix::InnerIterator it(A, outer); it; ++it) {
            triplets.emplace_back(it.row(),     it.col(),     it.value());
            triplets.emplace_back(it.row() + n, it.col() + n, it.value());
        }
    }

    for (Index outer = 0; outer < static_cast<Index>(B.outerSize()); ++outer) {
        for (SparseMatrix::InnerIterator it(B, outer); it; ++it) {
            triplets.emplace_back(it.row(),     it.col() + n, -it.value());
            triplets.emplace_back(it.row() + n, it.col(),      it.value());
        }
    }

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
 * Solves the harmonic response for each configured excitation frequency.
 *
 * Supports and topology constraints are assembled from the current ModelData.
 * The constrained stiffness/mass operators and Rayleigh damping define the real
 * block representation of (K - omega^2 M + i omega C) u_hat = f_hat. The model
 * evaluates current direct loads and collector-activated conditions at each frequency,
 * preserving existing amplitude semantics and generalized-force transformations.
 *
 * Complex displacement components are expanded into global nodal fields for
 * stress, strain and requested harmonic output. This analysis owns the frequency
 * sweep and solver parameters, while condition history remains unchanged.
 */
void LinearHarmonic::run() {
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR HARMONIC RESPONSE ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

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

    model->assign_sections();
    model->step_begin();

    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_structural_dof_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat); },
        "building constraints");

    report_constraint_groups(groups);
    auto equations = groups.flatten();

    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K");

    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M");

    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                K.rows(),
                options);
        },
        "building constraint transformer");

    logging::error(transformer->homogeneous(),
        "LinearHarmonic: harmonic constraints must be homogeneous");
    logging::error(transformer->feasible(),
        "LinearHarmonic: constraint system is infeasible");

    auto Kr = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(K); },
        "assembling reduced stiffness matrix Kr");

    auto Mr = Timer::measure(
        [&]() { return transformer->reduce_secondary_matrix(M); },
        "assembling reduced mass matrix Mr");

    auto Cr = Timer::measure(
        [&]() { return damping.build(Mr, Kr); },
        "constructing reduced Rayleigh damping matrix Cr");

    DynamicMatrix basis(Kr.rows(), 0);
    DynamicVector eigenvalues;
    DynamicMatrix modal_stiffness;
    DynamicMatrix modal_mass;
    DynamicMatrix modal_damping;
    if (modal_basis) {
        const Precision upper_frequency  = 2 * *std::max_element(frequencies.begin(), frequencies.end());
        const Precision upper_omega      = 2 * pi * upper_frequency;
        const Precision upper_eigenvalue = upper_omega * upper_omega;
        logging::error(std::isfinite(upper_eigenvalue),
            "LinearHarmonic: modal basis upper eigenvalue must be finite");

        if (upper_eigenvalue > 0 && Kr.rows() > 1) {
            solver::EigvalOpts opts;
            opts.mode = solver::EigvalMode::ShiftInvert;
            opts.sort = solver::EigvalOpts::Sort::LargestMagn;
            auto pairs = Timer::measure(
                [&]() { return solver::eigvals(device, Kr, Mr, Precision(0), upper_eigenvalue, opts); },
                "constructing harmonic modal basis");

            const Index num_modes = static_cast<Index>(pairs.size());
            basis.resize(Kr.rows(), num_modes);
            eigenvalues.resize(num_modes);
            for (Index mode = 0; mode < num_modes; ++mode) {
                const auto& pair = pairs[static_cast<std::size_t>(mode)];
                const Precision modal_mass = pair.vector.dot(Mr * pair.vector);
                logging::error(std::isfinite(modal_mass) && modal_mass > 0,
                    "LinearHarmonic: modal mass must be finite and positive");
                basis.col(mode)   = pair.vector / std::sqrt(modal_mass);
                eigenvalues(mode) = pair.value;
            }
        }

        if (basis.cols() > 0) {
            auto load_basis = model->build_load_basis(loads);
            if (!load_basis.empty()) {
                const Index num_patterns = static_cast<Index>(load_basis.size());
                DynamicMatrix load_vectors(Kr.rows(), num_patterns);
                for (Index pattern = 0; pattern < num_patterns; ++pattern) {
                    const auto& field = load_basis[pattern].second;
                    const DynamicVector f = mattools::reduce_mat_to_vec(active_dof_idx_mat, field);
                    load_vectors.col(pattern) = transformer->assemble_system_rhs(K, f);
                }

                const DynamicMatrix static_responses = Timer::measure(
                    [&]() { return solver::solve(device, method, Kr, load_vectors, solver::DirectSolverMatrixType::SPD); },
                    "solving harmonic static load basis");
                logging::error(static_responses.allFinite(),
                    "LinearHarmonic: static load basis contains invalid values");

                const Precision residual_tolerance = Precision(1e-10);
                for (Index pattern = 0; pattern < num_patterns; ++pattern) {
                    DynamicVector residual = static_responses.col(pattern);
                    const Precision static_mass = residual.dot(Mr * residual);
                    logging::error(std::isfinite(static_mass) && static_mass >= 0,
                        "LinearHarmonic: static response mass norm must be finite and non-negative");
                    if (static_mass == 0) {
                        continue;
                    }

                    for (int pass = 0; pass < 2; ++pass) {
                        const DynamicVector projection = basis.transpose() * (Mr * residual);
                        residual -= basis * projection;
                    }

                    const Precision residual_mass = residual.dot(Mr * residual);
                    logging::error(std::isfinite(residual_mass) && residual_mass >= 0,
                        "LinearHarmonic: residual mass norm must be finite and non-negative");
                    if (residual_mass <= residual_tolerance * residual_tolerance * static_mass) {
                        continue;
                    }

                    const Eigen::Index column = basis.cols();
                    basis.conservativeResize(Kr.rows(), column + 1);
                    basis.col(column) = residual / std::sqrt(residual_mass);
                }
            }

            modal_stiffness = basis.transpose() * (Kr * basis);
            modal_mass     = basis.transpose() * (Mr * basis);
            modal_damping  = damping.alpha * modal_mass + damping.beta * modal_stiffness;
        }

        logging::info(true             , "Modal basis modes          : ", eigenvalues.size());
        logging::info(true             , "Residual basis vectors     : ", basis.cols() - eigenvalues.size());
        logging::info(true             , "Modal basis upper frequency: ", upper_frequency);
        logging::info(basis.cols() == 0, "Modal basis is empty; using the direct harmonic solver");
    }

    writer->add_loadcase(id, io::writer::WriterStepType::Dynamic);

    logging::info(true, "Frequency sweep");
    logging::info(true,
                  std::setw(6),  "Idx",
                  std::setw(16), "Frequency",
                  std::setw(18), "Response norm");
    logging::info(true, "----------------------------------------");

    Timer sweep_timer;
    sweep_timer.start();
    for (Index i = 0; i < static_cast<Index>(frequencies.size()); ++i) {
        const Precision frequency = frequencies[i];
        const Precision omega     = 2.0 * pi * frequency;

        DynamicVector u_real;
        DynamicVector u_imag;

        run_quiet([&]() {
            // Abaqus-style amplitudes use frequency as STEP TIME in a frequency
            // domain procedure. Rebuild only the load vector for each frequency;
            // K, M, C and the constraint transformation remain unchanged.
            model::Field global_load_mat = model->build_load_matrix(frequency);
            const DynamicVector f = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat);
            const DynamicVector fr = transformer->assemble_system_rhs(K, f);

            DynamicVector q_real;
            DynamicVector q_imag;
            if (basis.cols() > 0) {
                const Index num_modes = static_cast<Index>(basis.cols());
                const DynamicMatrix A = modal_stiffness - omega * omega * modal_mass;
                const DynamicMatrix B = omega * modal_damping;

                DynamicMatrix block_matrix(2 * num_modes, 2 * num_modes);
                block_matrix.topLeftCorner    (num_modes, num_modes) =  A;
                block_matrix.topRightCorner   (num_modes, num_modes) = -B;
                block_matrix.bottomLeftCorner (num_modes, num_modes) =  B;
                block_matrix.bottomRightCorner(num_modes, num_modes) =  A;

                DynamicVector rhs = DynamicVector::Zero(2 * num_modes);
                rhs.head(num_modes) = basis.transpose() * fr;

                Eigen::FullPivLU<DynamicMatrix> modal_solver(block_matrix);
                logging::error(modal_solver.isInvertible(),
                    "LinearHarmonic: projected dynamic stiffness is singular");
                const DynamicVector solution = modal_solver.solve(rhs);
                logging::error(solution.allFinite(),
                    "LinearHarmonic: projected harmonic solution contains invalid values");
                q_real = basis * solution.head(num_modes);
                q_imag = basis * solution.tail(num_modes);
            } else {
                SparseMatrix A = Kr - omega * omega * Mr;
                SparseMatrix B = omega * Cr;
                A.makeCompressed();
                B.makeCompressed();

                SparseMatrix block_matrix = build_block_matrix(A, B);
                const DynamicVector zero = DynamicVector::Zero(fr.size());
                const DynamicVector rhs  = build_block_rhs(fr, zero);

                DynamicVector solution = solve(
                    device,
                    method,
                    block_matrix,
                    rhs,
                    solver::DirectSolverMatrixType::General);

                const Index n = static_cast<Index>(fr.size());
                q_real = solution.head(n);
                q_imag = solution.tail(n);
            }

            u_real = transformer->recover_displacement(q_real);
            u_imag = transformer->recover_displacement(q_imag);

            auto displacement_real = mattools::expand_vec_to_mat(active_dof_idx_mat, u_real);
            auto displacement_imag = mattools::expand_vec_to_mat(active_dof_idx_mat, u_imag);

            // Harmonic displacement is solved directly as real and imaginary
            // components. Stress/strain components are recovered lazily from
            // the corresponding displacement component only when requested.
            using io::writer::OutputField;
            output.begin_frame(frequency, "_" + std::to_string(i + 1));
            output.provide(OutputField::DISPLACEMENT_REAL,    displacement_real);
            output.provide(OutputField::DISPLACEMENT_IMAG,    displacement_imag);
            output.provide(OutputField::EXTERNAL_FORCES,      global_load_mat);
            output.write_frame(*writer, model->_data.get());
        });

        logging::info(true,
                      std::setw(6),  i + 1,
                      std::setw(16), std::fixed, std::setprecision(6), frequency,
                      std::setw(18), std::scientific, std::setprecision(6),
                      std::sqrt(u_real.squaredNorm() + u_imag.squaredNorm()));
    }
    sweep_timer.stop();
    logging::info(true, "Harmonic frequency sweep elapsed: ", sweep_timer.elapsed(), " ms");

    logging::info(true, "");
    model->step_end();
}

} // namespace loadcase
} // namespace fem
