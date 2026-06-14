/**
 * @file nonlinear_static.cpp
 * @brief Implements incremental nonlinear static analysis.
 */

#include "nonlinear_static.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../core/config.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"
#include "../writer/write_mtx.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

void set_positions(model::Model& model, const model::Field& positions) {
    logging::error(model._data != nullptr, "NonlinearStatic: model data not initialized");
    logging::error(model._data->positions != nullptr, "NonlinearStatic: positions field not initialized");
    *model._data->positions = positions;
}

model::Field current_positions_from_displacement(const model::Field& reference,
                                                 const model::Field& displacement) {
    logging::error(reference.domain == model::FieldDomain::NODE,
                   "NonlinearStatic: reference positions must use NODE domain");
    logging::error(displacement.domain == model::FieldDomain::NODE,
                   "NonlinearStatic: displacement must use NODE domain");
    logging::error(reference.rows == displacement.rows,
                   "NonlinearStatic: reference/displacement row mismatch");

    model::Field current = reference;
    for (Index i = 0; i < current.rows; ++i) {
        for (Index d = 0; d < std::min<Index>(3, current.components); ++d) {
            current(i, d) = reference(i, d) + displacement(i, d);
        }
    }
    return current;
}

model::Field subtract_field(const model::Field& lhs, const model::Field& rhs, const std::string& name) {
    model::Field result = lhs;
    result.name = name;
    result -= rhs;
    return result;
}

DynamicVector homogeneous_increment(const ConstraintTransformer& transformer,
                                    const DynamicVector& dq) {
    return transformer.map().T() * dq;
}

void check_finite(const SparseMatrix& matrix, const char* name) {
    bool bad = false;
    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            if (!std::isfinite(it.value())) {
                bad = true;
                break;
            }
        }
        if (bad) break;
    }
    logging::error(!bad, name, " contains NaN/Inf entries");
}

int apply_zero_stiffness_regularization(SparseMatrix& matrix, Precision alpha) {
    if (alpha <= Precision(0) || matrix.rows() == 0) {
        return 0;
    }

    Precision diagonal_sum = Precision(0);
    Index diagonal_count = 0;
    Precision entry_sum = Precision(0);
    Index entry_count = 0;
    DynamicVector row_norm = DynamicVector::Zero(matrix.rows());

    for (int col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(matrix, col); it; ++it) {
            const Precision value_abs = std::abs(it.value());
            row_norm(it.row()) += value_abs;
            if (value_abs > Precision(0)) {
                entry_sum += value_abs;
                ++entry_count;
            }
            if (it.row() == it.col() && value_abs > Precision(0)) {
                diagonal_sum += value_abs;
                ++diagonal_count;
            }
        }
    }

    Precision stiffness_scale = Precision(0);
    if (diagonal_count > 0) {
        stiffness_scale = diagonal_sum / static_cast<Precision>(diagonal_count);
    } else if (entry_count > 0) {
        stiffness_scale = entry_sum / static_cast<Precision>(entry_count);
    }

    if (stiffness_scale <= Precision(0) || !std::isfinite(stiffness_scale)) {
        return 0;
    }

    const Precision min_row_stiffness = std::max(
        alpha * stiffness_scale,
        std::numeric_limits<Precision>::epsilon());

    TripletList additions;
    for (int row = 0; row < matrix.rows(); ++row) {
        if (row_norm(row) < min_row_stiffness) {
            additions.emplace_back(row, row, min_row_stiffness - row_norm(row));
        }
    }

    if (additions.empty()) {
        return 0;
    }

    SparseMatrix regularization(matrix.rows(), matrix.cols());
    regularization.setFromTriplets(additions.begin(), additions.end());
    matrix += regularization;
    matrix.makeCompressed();
    return static_cast<int>(additions.size());
}

std::string nonlinear_solver_backend(solver::SolverDevice device) {
    std::ostringstream out;
    if (device == solver::GPU) {
#ifdef SUPPORT_GPU
#ifdef USE_CUDSS
        out << "GPU DIRECT cuDSS";
#else
        out << "GPU DIRECT CUDA fallback";
#endif
#else
        out << "CPU DIRECT";
#ifdef USE_MKL
        out << " MKL PardisoLDLT";
#else
        out << " Eigen SimplicialLDLT";
#endif
        out << " (GPU requested, build has no GPU support)";
#endif
    } else {
        out << "CPU DIRECT";
#ifdef USE_MKL
        out << " MKL PardisoLDLT";
#else
        out << " Eigen SimplicialLDLT";
#endif
    }
    out << ", max threads = " << global_config.max_threads;
    return out.str();
}

DynamicVector solve_quietly(solver::SolverDevice device,
                            solver::SolverMethod method,
                            SparseMatrix& matrix,
                            DynamicVector& rhs) {
    const bool logging_was_enabled = logging::is_enabled();
    if (logging_was_enabled) {
        logging::disable();
    }

    try {
        DynamicVector result = solver::solve(device, method, matrix, rhs);
        if (logging_was_enabled) {
            logging::enable();
        }
        return result;
    } catch (...) {
        if (logging_was_enabled) {
            logging::enable();
        }
        throw;
    }
}

} // namespace

NonlinearStatic::NonlinearStatic(ID id, reader::Writer* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

void NonlinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "NONLINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    logging::error(num_increments > 0, "NONLINEARSTATIC requires INCREMENTS > 0");
    logging::error(max_iterations > 0, "NONLINEARSTATIC requires MAXITER > 0");
    logging::error(tolerance > Precision(0), "NONLINEARSTATIC requires TOL > 0");
    logging::error(constraint_method == ConstraintTransformer::Method::NullSpace,
                   "NONLINEARSTATIC currently supports only NULLSPACE constraints");
    logging::error(method == solver::DIRECT,
                   "NONLINEARSTATIC currently supports only DIRECT solver method");

    model->assign_sections();
    if (!model->_data->element_ip_offsets) {
        model->_data->initialize_element_enumeration();
    }

    model::Field reference_positions = *model->_data->positions;

    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto global_load_total = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing total load matrix (node x 6)");

    auto f_total = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_total); },
        "reducing total load matrix -> active RHS vector f");

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints");
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    const Index n_active = active_dof_idx_mat.maxCoeff() + 1;
    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions options;
            options.method = constraint_method;
            options.set.scale_columns = true;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                n_active,
                options);
        },
        "building constraint transformer");

    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)     : ", transformer->report().m);
    logging::info(true, "n (cols of C)     : ", transformer->report().n);
    logging::info(true, "rank(C)           : ", transformer->rank());
    logging::info(true, "method            : ", transformer->method_name());
    logging::info(true, "masters (n-r)     : ", transformer->n_master());
    logging::info(true, "homogeneous       : ", transformer->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", transformer->feasible() ? "true" : "false");
    if (!transformer->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", transformer->report().residual_norm);
    }
    logging::down();

    DynamicVector q_zero = DynamicVector::Zero(transformer->n_master());
    DynamicVector u_total = transformer->recover_u(q_zero);

    model::Field displacement =
        mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);
    model::Field current_positions =
        current_positions_from_displacement(reference_positions, displacement);

    model::Field final_ip_stress{"IP_STRESS", model::FieldDomain::ELEMENT_IP, 1, 6};
    final_ip_stress.set_zero();
    model::Field final_internal{"INTERNAL_FORCES", model::FieldDomain::NODE,
                                static_cast<Index>(model->_data->max_nodes), 6};
    final_internal.set_zero();
    SparseMatrix final_Kt;
    SparseMatrix final_A;

    logging::info(true, "");
    logging::info(true, "Solver backend: ", nonlinear_solver_backend(device));
    logging::info(true, "");
    logging::info(true, " inc iter      lambda        rel_res          du_norm   asm_ms solve_ms");
    logging::info(true, "--------------------------------------------------------------------------");

    for (int increment = 1; increment <= num_increments; ++increment) {
        const Precision load_factor =
            static_cast<Precision>(increment) / static_cast<Precision>(num_increments);
        const DynamicVector f_ext = load_factor * f_total;

        bool converged = false;
        for (int iter = 1; iter <= max_iterations; ++iter) {
            Timer asm_timer;
            asm_timer.start();

            displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);

            set_positions(*model, reference_positions);
            auto ip_stress = Timer::measure(
                [&]() { return model->compute_stress_state(displacement, true); },
                "computing nonlinear stress state",
                false);

            current_positions = current_positions_from_displacement(reference_positions, displacement);
            set_positions(*model, current_positions);

            auto K = Timer::measure(
                [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
                "assembling current material stiffness K",
                false);
            auto Kg = Timer::measure(
                [&]() { return model->build_geom_stiffness_matrix(active_dof_idx_mat, ip_stress); },
                "assembling current geometric stiffness K_g",
                false);
            SparseMatrix Kt = K + Kg;
            Kt.makeCompressed();
            if (regularize_zero_stiffness_rows) {
                apply_zero_stiffness_regularization(Kt, zero_stiffness_regularization_alpha);
            }
            check_finite(Kt, "K_t");

            auto internal_mat = Timer::measure(
                [&]() { return model->build_internal_force_nonlinear(ip_stress); },
                "assembling nonlinear internal force",
                false);
            auto f_int = Timer::measure(
                [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, internal_mat); },
                "reducing internal force matrix -> active vector",
                false);

            DynamicVector residual = f_ext - f_int;
            DynamicVector reduced_residual;
            transformer->apply_Tt(residual, reduced_residual);

            DynamicVector reduced_external;
            transformer->apply_Tt(f_ext, reduced_external);
            const Precision denom = std::max<Precision>(reduced_external.norm(), Precision(1));
            const Precision rel_residual = reduced_residual.norm() / denom;
            Precision du_norm = Precision(0);
            Time solve_ms = Time(0);

            final_ip_stress = ip_stress;
            final_internal = internal_mat;
            final_Kt = Kt;

            if (rel_residual <= tolerance) {
                final_A = transformer->assemble_A(Kt);
                asm_timer.stop();
                logging::info(true,
                              std::setw(4), increment,
                              std::setw(5), iter,
                              std::scientific, std::setprecision(3),
                              std::setw(12), load_factor,
                              std::setw(15), rel_residual,
                              std::setw(15), du_norm,
                              std::fixed, std::setprecision(1),
                              std::setw(9), asm_timer.elapsed(),
                              std::setw(9), solve_ms);
                converged = true;
                break;
            }

            auto A = Timer::measure(
                [&]() { return transformer->assemble_A(Kt); },
                "assembling reduced tangent A = T^T K_t T",
                false);
            check_finite(A, "A");
            logging::error(reduced_residual.allFinite(), "Reduced residual contains NaN/Inf entries");
            asm_timer.stop();

            Timer solve_timer;
            solve_timer.start();
            auto dq = solve_quietly(device, method, A, reduced_residual);
            solve_timer.stop();
            solve_ms = solve_timer.elapsed();

            DynamicVector du = homogeneous_increment(*transformer, dq);
            du_norm = du.norm();
            u_total += du;

            final_A = A;

            logging::info(true,
                          std::setw(4), increment,
                          std::setw(5), iter,
                          std::scientific, std::setprecision(3),
                          std::setw(12), load_factor,
                          std::setw(15), rel_residual,
                          std::setw(15), du_norm,
                          std::fixed, std::setprecision(1),
                          std::setw(9), asm_timer.elapsed(),
                          std::setw(9), solve_ms);
        }

        logging::warning(converged,
                         "NONLINEARSTATIC increment ", increment,
                         " did not converge after ", max_iterations,
                         " iterations; continuing with current state");
    }

    displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);

    set_positions(*model, reference_positions);
    final_ip_stress = Timer::measure(
        [&]() { return model->compute_stress_state(displacement, true); },
        "computing final nonlinear stress state");

    auto [final_stress, final_strain] = Timer::measure(
        [&]() { return model->compute_stress_nodal(displacement, true); },
        "computing final nonlinear nodal stress/strain");

    auto [final_stress_top, final_stress_bot] = Timer::measure(
        [&]() { return model->compute_stress_top_bot(displacement, true); },
        "computing final nonlinear top/bottom stress");

    current_positions = current_positions_from_displacement(reference_positions, displacement);
    set_positions(*model, current_positions);

    auto final_K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "assembling final current material stiffness K");
    auto final_Kg = Timer::measure(
        [&]() { return model->build_geom_stiffness_matrix(active_dof_idx_mat, final_ip_stress); },
        "assembling final current geometric stiffness K_g");
    final_Kt = final_K + final_Kg;
    final_Kt.makeCompressed();
    if (regularize_zero_stiffness_rows) {
        apply_zero_stiffness_regularization(final_Kt, zero_stiffness_regularization_alpha);
    }
    final_A = transformer->assemble_A(final_Kt);
    final_internal = Timer::measure(
        [&]() { return model->build_internal_force_nonlinear(final_ip_stress); },
        "assembling final nonlinear internal force");

    auto global_load_final = global_load_total;
    auto reaction_full = subtract_field(final_internal, global_load_final, "REACTION_FORCES_RAW");

    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);
    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 && e.node_id < support_mask.rows() &&
                e.dof >= 0 && e.dof < support_mask.cols() &&
                active_dof_idx_mat(e.node_id, e.dof) != -1) {
                support_mask(e.node_id, e.dof) = true;
            }
        }
    }

    model::Field reaction_masked{"REACTION_FORCES", model::FieldDomain::NODE,
                                 reaction_full.rows,
                                 reaction_full.components};
    reaction_masked.fill_nan();
    for (int i = 0; i < reaction_masked.rows; ++i) {
        for (int j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = reaction_full(i, j);
            }
        }
    }

    writer->add_loadcase(id);
    writer->write_field(displacement     , "DISPLACEMENT", model->_data.get());
    writer->write_field(final_strain     , "STRAIN", model->_data.get());
    writer->write_field(final_stress     , "STRESS", model->_data.get());
    writer->write_field(final_stress_top , "STRESS_TOP", model->_data.get());
    writer->write_field(final_stress_bot , "STRESS_BOT", model->_data.get());
    writer->write_field(global_load_final, "EXTERNAL_FORCES", model->_data.get());
    writer->write_field(final_internal   , "INTERNAL_FORCES", model->_data.get());
    writer->write_field(reaction_masked  , "REACTION_FORCES", model->_data.get());

    if (!stiffness_file.empty()) {
        write_mtx(stiffness_file + "_Kt.mtx", final_Kt);
        write_mtx(stiffness_file + "_A.mtx", final_A);
    }

    set_positions(*model, reference_positions);
}

} // namespace loadcase
} // namespace fem
