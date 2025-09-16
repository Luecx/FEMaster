/******************************************************************************
 * @file LinearStaticTopo.cpp
 * @brief Linear static analysis with topology optimization (density/orientation).
 *
 * Pipeline (null-space constrained, like LinearStatic):
 *  - Assign sections; push SIMP-like element scaling (density^p) and orientation.
 *  - Build unconstrained DOF indexing (node×6 → active dof id or -1).
 *  - Build constraint equations (supports/ties/couplings → C u = d).
 *  - Assemble active stiffness K(ρ,p,θ) and global loads; reduce loads to f.
 *  - Build ConstraintTransformer → u = u_p + T q.
 *  - Reduce: A = Tᵀ K T, b = Tᵀ (f − K u_p); solve A q = b.
 *  - Recover u, reactions r = K u − f; expand to node×6.
 *  - Post: stress/strain, compliance & sensitivities, volumes; write results.
 *
 * @date 27.08.2024
 * @author Finn Eggers
 ******************************************************************************/

#include "linear_static_topo.h"

#include "../core/logging.h"
#include "../solve/eigen.h"

#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"
#include "../constraints/constraint_set.h"
#include "../constraints/constraint_builder.h"
#include "../constraints/constraint_map.h"

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"

#include <limits>

using fem::constraint::ConstraintTransformer;

namespace fem { namespace loadcase {

LinearStaticTopo::LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model)
    : LinearStatic(id, writer, model),
      density(model->_data->max_elems, 1),
      orientation(model->_data->max_elems, 3)
{
    density.setOnes();
    orientation.setZero();
}

void LinearStaticTopo::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC TOPO");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    m_model->assign_sections();

    // SIMP-like element scaling and per-element orientation for K assembly.
    const ElementData stiffness_scalar = density.array().pow(exponent);
    m_model->_data->create_data(model::ElementDataEntries::TOPO_STIFFNESS, 1);
    m_model->_data->create_data(model::ElementDataEntries::TOPO_ANGLES   , 3);
    m_model->_data->get(model::ElementDataEntries::TOPO_STIFFNESS) = stiffness_scalar;
    m_model->_data->get(model::ElementDataEntries::TOPO_ANGLES)    = orientation;

    // (1) Unconstrained DOF indexing
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraint equations
    auto equations = Timer::measure(
        [&]() { return this->m_model->build_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // (3) Global loads (node×6) and reduction to active RHS
    auto global_load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(loads); },
        "constructing load matrix (node×6)"
    );
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix → active RHS vector f"
    );

    // (4) Active stiffness with topo scaling/orientation
    auto K = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat, stiffness_scalar); },
        "constructing stiffness matrix K(ρ^p, θ)"
    );

    // (5) Constraint transformer
    ConstraintTransformer::BuildOptions copt;
    copt.set.scale_columns = true;
    ConstraintTransformer CT(
        equations,
        active_dof_idx_mat,
        K.rows(),
        copt
    );

    // (6) Reduced system A q = b
    auto A = Timer::measure([&]() { return CT.assemble_A(K); }, "assembling A = T^T K T");
    auto b = Timer::measure([&]() { return CT.assemble_b(K, f); }, "assembling b = T^T (f - K u_p)");

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b); },
        "solving reduced system A q = b"
    );

    // (7) Recover u, reactions r
    auto u = Timer::measure([&]() { return CT.recover_u(q); }, "recovering full displacement vector u");
    auto r = Timer::measure([&]() { return CT.reactions(K, f, q); }, "computing reactions r = K u - f");

    // (8) Expand to node×6
    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form"
    );
    auto global_force_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r); },
        "expanding reactions to matrix form"
    );

    // (9) Post-processing: stress/strain
    NodeData stress, strain;
    std::tie(stress, strain) = Timer::measure(
        [&]() { return m_model->compute_stress_strain(global_disp_mat); },
        "Interpolating stress and strain at nodes"
    );

    // (10) Topology metrics: compliance & gradients, volumes
    ElementData compliance_raw = m_model->compute_compliance(global_disp_mat);
    ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);
    ElementData volumes        = m_model->compute_volumes();
    ElementData angle_grad     = m_model->compute_compliance_angle_derivative(global_disp_mat);

    // (11) Write results
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(global_disp_mat , "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain          , "STRAIN");
    m_writer->write_eigen_matrix(stress          , "STRESS");
    m_writer->write_eigen_matrix(global_load_mat , "DOF_LOADS");
    m_writer->write_eigen_matrix(global_force_mat, "NODAL_FORCES");
    m_writer->write_eigen_matrix(compliance_raw  , "COMPLIANCE_RAW");
    m_writer->write_eigen_matrix(compliance_adj  , "COMPLIANCE_ADJ");
    m_writer->write_eigen_matrix(dens_grad       , "DENS_GRAD");
    m_writer->write_eigen_matrix(volumes         , "VOLUME");
    m_writer->write_eigen_matrix(density         , "DENSITY");
    m_writer->write_eigen_matrix(angle_grad      , "ORIENTATION_GRAD");
    m_writer->write_eigen_matrix(orientation     , "ORIENTATION");

    // (12) Cleanup element scratch data
    m_model->_data->remove(model::ElementDataEntries::TOPO_STIFFNESS);
    m_model->_data->remove(model::ElementDataEntries::TOPO_ANGLES);

    // (13) Diagnostics (optional)
    {
        DynamicVector resid = K * u - f;
        DynamicVector red   = CT.map().apply_Tt(resid);
        logging::info(true, "");
        logging::info(true, "Post-checks");
        logging::up();
        logging::info(true, "||u||2                : ", u.norm());
        logging::info(true, "||C u - d||2          : ", (CT.set().C * u - CT.set().d).norm());
        logging::info(true, "||K u - f||2          : ", resid.norm());
        logging::info(true, "||T^T (K u - f)||2    : ", red.norm());
        logging::down();
    }
}

}} // namespace fem::loadcase
