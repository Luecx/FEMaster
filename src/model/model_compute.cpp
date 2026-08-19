/**
 * @file model_compute.cpp
 * @brief Implements model-wide finite-element result recovery.
 *
 * The routines allocate result fields in their physical storage domains and
 * delegate element-specific evaluation to compiled structural elements.
 * Integration-point results use the offsets established during compilation;
 * element-nodal quantities are recovered into global nodal fields through the
 * weighting operation provided by `ModelData`.
 *
 * Nodal stress recovery, shell-face recovery, beam section forces and shear
 * flow may evaluate elements with OpenMP. Eigen and MKL internal threading are
 * temporarily restricted during those loops to avoid nested oversubscription.
 *
 * @see Model
 * @see ModelData
 * @see StructuralElement
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../core/config.h"
#include "element/element_structural.h"
#include "model.h"
#include "shell/s8.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem::model {

/**
 * Recovers the complete structural stress state at all integration points.
 *
 * The compiled element-IP offset field maps every structural element to its
 * contiguous range in the returned `ELEMENT_IP` field. Each element evaluates
 * its formulation-specific stress state from the supplied displacement field
 * and the requested linear or Green-Lagrange nonlinear kinematics.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @param use_green_lagrange_nl Enables the nonlinear Green-Lagrange strain path
 *                              in supporting element formulations.
 * @return Integration-point stress-state field with the established
 *         eight-component FEMaster layout.
 */
Field Model::compute_stress_state(Field& displacement, bool use_green_lagrange_nl) {
    // Validate and access the compiled integration-point enumeration
    logging::error(_data->element_ip_offsets != nullptr,
        "element IP offset field has not been initialized");

    const auto& ip_enum = *_data->element_ip_offsets;
    const Index element_count = static_cast<Index>(_data->elements.size());
    const Index total_ips = _data->field_rows(FieldDomain::ELEMENT_IP);

    // Allocate one zeroed row for every compiled element integration point
    Field ip_stress{"IP_STRESS", FieldDomain::ELEMENT_IP, total_ips, 8};
    ip_stress.set_zero();

    // Delegate stress-state recovery to each compiled structural element
    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && static_cast<Index>(eid) < element_count,
                "Element id out of range in compute_stress_state: ", eid);

            const Index ip_offset = static_cast<Index>(ip_enum(static_cast<Index>(eid), 0));
            logging::error(ip_offset <= total_ips,
                "Invalid IP offset for element ", eid, ": ", ip_offset,
                " / total=", total_ips);

            sel->compute_stress_state(
                ip_stress,
                displacement,
                static_cast<int>(ip_offset),
                use_green_lagrange_nl
            );
        }
    }

    // Reject invalid constitutive or kinematic results before returning the field
    ip_stress.check_finite("Stress state");
    return ip_stress;
}

/**
 * Recovers averaged nodal stress and strain fields.
 *
 * Every supporting structural element supplies its natural nodal recovery
 * coordinates and writes six-component stress and strain values into
 * `ELEMENT_NODAL` storage. The compiled element-nodal offsets keep parallel
 * writes disjoint. `ModelData::element_nodal_to_nodal()` then averages the
 * contributing element values at each global node using unit element weights.
 *
 * Elements without nodal recovery coordinates do not participate.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @param use_green_lagrange_nl Enables Green-Lagrange strain recovery in
 *                              supporting nonlinear formulations.
 * @return Pair containing the global nodal stress field followed by strain.
 */
std::tuple<Field, Field> Model::compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate disjoint element-nodal recovery fields and participation weights
    Field element_stress {"ELEMENT_NODAL_STRESS", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_strain {"ELEMENT_NODAL_STRAIN", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_ELEMENT_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    element_stress.set_zero();
    element_strain.set_zero();
    element_weights.set_zero();

    // Prevent nested Eigen/MKL threading inside the element-parallel recovery loop
    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    // Recover element-nodal values into non-overlapping compiled row ranges
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            RowMatrix rst = sel->stress_strain_nodal_rst();
            if (rst.rows() == 0) continue;

            logging::error(rst.rows() == sel->n_nodes(),
                "Element ", sel->elem_id, " returned ", rst.rows(),
                " nodal stress coordinates, expected ", sel->n_nodes());
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_stress_strain(
                &element_strain,
                &element_stress,
                displacement,
                rst,
                static_cast<int>(offset),
                use_green_lagrange_nl
            );
            element_weights(static_cast<Index>(sel->elem_id), 0) = Precision(1);
        }
    }

    // Restore the configured linear-algebra thread counts after element recovery
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);

    // Average participating element values onto global nodes and validate them
    Field stress = _data->element_nodal_to_nodal(element_stress, element_weights, "STRESS");
    Field strain = _data->element_nodal_to_nodal(element_strain, element_weights, "STRAIN");
    stress.check_finite("Nodal stress");
    strain.check_finite("Nodal strain");
    return {stress, strain};
}

/**
 * Recovers averaged stress on the top and bottom shell faces.
 *
 * Shell recovery coordinates use the formulation-provided natural nodal points
 * with the through-thickness coordinate set to `+1` for the top face and `-1`
 * for the bottom face. Non-shell structural elements use their unchanged nodal
 * recovery coordinates for both fields. Element-nodal results are averaged onto
 * global nodes using the same participation weighting as regular nodal stress.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @param use_green_lagrange_nl Enables Green-Lagrange strain kinematics while
 *                              evaluating the requested stress values.
 * @return Pair containing global nodal top-face and bottom-face stress fields.
 */
std::tuple<Field, Field> Model::compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate separate element-nodal fields for both faces and their weights
    Field element_top    {"ELEMENT_NODAL_STRESS_TOP", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_bot    {"ELEMENT_NODAL_STRESS_BOT", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_TOP_BOT_ELEMENT_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    element_top.set_zero();
    element_bot.set_zero();
    element_weights.set_zero();

    // Prevent nested Eigen/MKL threading inside the element-parallel recovery loop
    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    // Evaluate both faces in each element's disjoint compiled row range
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            RowMatrix base_rst = sel->stress_strain_nodal_rst();
            if (base_rst.rows() == 0) continue;

            RowMatrix rst_bot = base_rst;
            RowMatrix rst_top = base_rst;

            // Select the two natural thickness faces for shell formulations
            if (sel->is_shell()) {
                for (int i = 0; i < base_rst.rows(); ++i) {
                    rst_bot(i, 2) = -1;
                    rst_top(i, 2) =  1;
                }
            }

            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_stress_strain(nullptr, &element_bot, displacement, rst_bot, static_cast<int>(offset), use_green_lagrange_nl);
            sel->compute_stress_strain(nullptr, &element_top, displacement, rst_top, static_cast<int>(offset), use_green_lagrange_nl);
            element_weights(static_cast<Index>(sel->elem_id), 0) = Precision(1);
        }
    }

    // Restore the configured linear-algebra thread counts after element recovery
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);

    // Average face values onto global nodes and reject invalid numerical output
    Field stress_top = _data->element_nodal_to_nodal(element_top, element_weights, "STRESS_TOP");
    Field stress_bot = _data->element_nodal_to_nodal(element_bot, element_weights, "STRESS_BOT");
    stress_top.check_finite("Nodal top stress");
    stress_bot.check_finite("Nodal bottom stress");
    return {stress_top, stress_bot};
}

/**
 * Recovers nodal shell force and moment resultants.
 *
 * Supporting elements accumulate their eight-component resultant convention and
 * one participation count at every connected global node. The model divides the
 * accumulated values by those counts to obtain an arithmetic nodal average;
 * nodes without shell contributions retain zero rows.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @return Averaged eight-component nodal shell-resultant field.
 */
Field Model::compute_shell_resultants(Field& displacement) {
    // Allocate nodal accumulators for resultants and contribution counts
    const Index node_count = _data->field_rows(FieldDomain::NODE);
    Field resultants{"SHELL_RESULTANTS", FieldDomain::NODE, node_count, 8};
    Field count{"SHELL_RESULTANTS_COUNT", FieldDomain::NODE, node_count, 1};
    resultants.set_zero();
    count.set_zero();

    // Accumulate formulation-specific shell resultants from structural elements
    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_shell_section_forces(resultants, count, displacement);
        }
    }

    // Convert accumulated contributions into arithmetic nodal averages
    for (Index i = 0; i < node_count; ++i) {
        if (count(i, 0) != Precision(0)) {
            for (Index j = 0; j < resultants.components; ++j) {
                resultants(i, j) /= count(i, 0);
            }
        }
    }

    // Validate the averaged result field before returning it
    resultants.check_finite("Shell resultants");
    return resultants;
}

/**
 * Computes one scalar compliance value for every compiled element.
 *
 * Structural formulations write their displacement-dependent compliance into
 * the row addressed by the dense element identifier. Non-structural or null
 * element slots retain zero.
 *
 * @param displacement Global nodal displacement field.
 * @return Scalar `ELEMENT`-domain compliance field.
 */
Field Model::compute_compliance(Field& displacement) {
    // Allocate the element-domain result with zeros for unsupported entities
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field compliance{"COMPLIANCE", FieldDomain::ELEMENT, element_count, 1};
    compliance.set_zero();

    // Delegate compliance evaluation to every structural element
    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_compliance(displacement, compliance);
        }
    }
    return compliance;
}

/**
 * Computes compliance derivatives with respect to material orientation angles.
 *
 * Every structural element writes its three orientation-gradient components to
 * its dense element row. Elements that do not implement an orientation
 * sensitivity leave the initialized zero row unchanged.
 *
 * @param displacement Global nodal displacement field.
 * @return Three-component `ELEMENT`-domain orientation-gradient field.
 */
Field Model::compute_compliance_angle_derivative(Field& displacement) {
    // Allocate the element-domain orientation gradient
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field results{"ORIENTATION_GRAD", FieldDomain::ELEMENT, element_count, 3};
    results.set_zero();

    // Delegate sensitivity evaluation to every structural element
    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_compliance_angle_derivative(displacement, results);
        }
    }
    return results;
}

/**
 * Computes the physical measure represented by every structural element.
 *
 * The formulation-specific `volume()` operation supplies the relevant scalar
 * element measure and stores it at the dense element identifier. Unsupported or
 * null entities retain zero.
 *
 * @return Scalar `ELEMENT`-domain volume field.
 */
Field Model::compute_volumes() {
    // Allocate the element-domain physical measure field
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field volumes{"VOLUME", FieldDomain::ELEMENT, element_count, 1};
    volumes.set_zero();

    // Evaluate the physical measure of every compiled structural element
    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            volumes(static_cast<Index>(el->elem_id)) = sel->volume();
        }
    }
    return volumes;
}

/**
 * Recovers six-component beam section forces at element nodes.
 *
 * Each structural formulation writes supported section-force quantities into
 * its disjoint range of the compiled `ELEMENT_NODAL` field. Element evaluation
 * may run in parallel while Eigen and MKL internal threading are restricted to
 * avoid nested oversubscription.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @return Six-component element-nodal beam section-force field.
 */
Field Model::compute_section_forces(Field& displacement) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate zeroed storage for every compiled element-node row
    Field beam_forces{"BEAM_SECTION_FORCES", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    beam_forces.set_zero();

    // Prevent nested Eigen/MKL threading inside the element-parallel recovery loop
    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    // Recover section forces into each element's disjoint nodal row range
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_beam_section_forces(beam_forces, displacement, static_cast<int>(offset));
        }
    }

    // Restore the configured linear-algebra thread counts
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);
    return beam_forces;
}

/**
 * Recovers scalar shear flow at element nodes.
 *
 * Supporting structural formulations write one shear-flow value into every
 * element-local nodal row addressed by the compiled offset field. Parallel
 * element evaluation uses the same nested-thread suppression as the other
 * element-nodal recovery operations.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @return Scalar `ELEMENT_NODAL` shear-flow field.
 */
Field Model::compute_shear_flow(Field& displacement) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate zeroed storage for every compiled element-node row
    Field shear_flow{"SHEAR_FLOW", FieldDomain::ELEMENT_NODAL, total_element_nodes, 1};
    shear_flow.set_zero();

    // Prevent nested Eigen/MKL threading inside the element-parallel recovery loop
    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    // Recover shear flow into each element's disjoint nodal row range
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_shear_flow(shear_flow, displacement, static_cast<int>(offset));
        }
    }

    // Restore the configured linear-algebra thread counts
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);
    return shear_flow;
}

}
