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
 * with A = K - omega^2 M and B = omega C. The first implementation assumes
 * real load amplitudes, hence f_i = 0, and Rayleigh damping C = alpha M + beta K.
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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

constexpr Precision pi = Precision(3.141592653589793238462643383279502884L);

SparseMatrix build_block_matrix(const SparseMatrix& A, const SparseMatrix& B) {
    logging::error(A.rows() == A.cols(), "LinearHarmonic: A must be square");
    logging::error(B.rows() == B.cols(), "LinearHarmonic: B must be square");
    logging::error(A.rows() == B.rows(), "LinearHarmonic: A and B size mismatch");

    using Triplet = Eigen::Triplet<Precision, Index>;

    const Index n = static_cast<Index>(A.rows());
    std::vector<Triplet> triplets;
    triplets.reserve(2 * (A.nonZeros() + B.nonZeros()));

    for (Index outer = 0; outer < A.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator it(A, outer); it; ++it) {
            triplets.emplace_back(it.row(),     it.col(),     it.value());
            triplets.emplace_back(it.row() + n, it.col() + n, it.value());
        }
    }

    for (Index outer = 0; outer < B.outerSize(); ++outer) {
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
    logging::error(real.size() == imag.size(), "LinearHarmonic: RHS size mismatch");

    DynamicVector result(2 * real.size());
    result.head(real.size()) = real;
    result.tail(real.size()) = imag;
    return result;
}

} // anonymous namespace

LinearHarmonic::LinearHarmonic(ID id, io::writer::ResultWriters* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

void LinearHarmonic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR HARMONIC RESPONSE ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    logging::error(!frequencies.empty(), "LinearHarmonic: no frequencies defined");
    logging::error(
        std::all_of(frequencies.begin(), frequencies.end(), [](Precision frequency) {
            return std::isfinite(frequency) && frequency >= 0.0;
        }),
        "LinearHarmonic: frequencies must be finite and non-negative");
    logging::error(
        constraint_method == ConstraintTransformer::Method::NullSpace,
        "LinearHarmonic: only NULLSPACE constraints are currently supported");

    model->assign_sections();
    model->step_begin();

    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing harmonic load-amplitude matrix");

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints");

    report_constraint_groups(groups);
    auto equations = groups.flatten();

    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K");

    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M");

    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing harmonic load matrix -> active RHS vector");

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

    auto fr = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(K, f); },
        "assembling reduced harmonic RHS");

    SparseMatrix Cr = rayleigh_alpha * Mr + rayleigh_beta * Kr;
    Cr.makeCompressed();

    writer->add_loadcase(id, io::writer::WriterStepType::Dynamic);

    logging::info(true, "");
    logging::info(true,
                  std::setw(8),  "Index",
                  std::setw(18), "Frequency",
                  std::setw(18), "Omega",
                  std::setw(18), "||u||");

    for (Index i = 0; i < static_cast<Index>(frequencies.size()); ++i) {
        const Precision frequency = frequencies[i];
        const Precision omega     = 2.0 * pi * frequency;

        SparseMatrix A = Kr - omega * omega * Mr;
        SparseMatrix B = omega * Cr;
        A.makeCompressed();
        B.makeCompressed();

        auto block_matrix = Timer::measure(
            [&]() { return build_block_matrix(A, B); },
            "constructing real harmonic block matrix");

        const DynamicVector zero = DynamicVector::Zero(fr.size());
        const DynamicVector rhs  = build_block_rhs(fr, zero);

        auto solution = Timer::measure(
            [&]() {
                return solve(device,
                             method,
                             block_matrix,
                             rhs,
                             solver::DirectSolverMatrixType::General);
            },
            "solving harmonic block system");

        const Index n = static_cast<Index>(fr.size());
        const DynamicVector q_real = solution.head(n);
        const DynamicVector q_imag = solution.tail(n);

        auto u_real = Timer::measure(
            [&]() { return transformer->recover_displacement(q_real); },
            "recovering real displacement vector");

        auto u_imag = Timer::measure(
            [&]() { return transformer->recover_displacement(q_imag); },
            "recovering imaginary displacement vector");

        auto displacement_real = Timer::measure(
            [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u_real); },
            "expanding real displacement vector");

        auto displacement_imag = Timer::measure(
            [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u_imag); },
            "expanding imaginary displacement vector");

        displacement_real.name = "DISPLACEMENT_REAL";
        displacement_imag.name = "DISPLACEMENT_IMAG";

        auto [stress_real, strain_real] = Timer::measure(
            [&]() { return model->compute_stress_nodal(displacement_real, false); },
            "computing real stress and strain");

        auto [stress_imag, strain_imag] = Timer::measure(
            [&]() { return model->compute_stress_nodal(displacement_imag, false); },
            "computing imaginary stress and strain");

        stress_real.name = "STRESS_REAL";
        stress_imag.name = "STRESS_IMAG";
        strain_real.name = "STRAIN_REAL";
        strain_imag.name = "STRAIN_IMAG";

        Timer::measure(
            [&]() {
                writer->write_field(displacement_real, displacement_real.name, model->_data.get(), frequency);
                writer->write_field(displacement_imag, displacement_imag.name, model->_data.get(), frequency);
                writer->write_field(stress_real, stress_real.name, model->_data.get(), frequency);
                writer->write_field(stress_imag, stress_imag.name, model->_data.get(), frequency);
                writer->write_field(strain_real, strain_real.name, model->_data.get(), frequency);
                writer->write_field(strain_imag, strain_imag.name, model->_data.get(), frequency);
            },
            "writing harmonic result fields");

        logging::info(true,
                      std::setw(8),  i + 1,
                      std::setw(18), std::fixed, std::setprecision(6), frequency,
                      std::setw(18), std::fixed, std::setprecision(6), omega,
                      std::setw(18), std::scientific, std::setprecision(6),
                      std::sqrt(u_real.squaredNorm() + u_imag.squaredNorm()));
    }

    model->step_end();
}

} // namespace loadcase
} // namespace fem
