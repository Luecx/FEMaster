/**
 * @file model_compute.cpp
 * @brief Implements model-wide finite-element result recovery.
 *
 * The routines allocate result fields in their physical storage domains and
 * delegate element-specific evaluation to compiled structural elements.
 * Integration-point results use the offsets established during compilation;
 * element-nodal quantities are projected to global nodes only when required,
 * using Field::to_nodal() to omit undefined element-local values.
 *
 * Nodal stress recovery, shell-face recovery, shell-resultant recovery, beam
 * section forces, shear flow and thermal heat-flux recovery may evaluate
 * elements with OpenMP. Eigen and MKL internal threading are temporarily
 * restricted during those loops to avoid nested oversubscription. Thermal flux
 * remains in disjoint ELEMENT_NODAL storage until nodal output requests
 * averaging through Field::to_nodal().
 *
 * @see Model
 * @see ModelData
 * @see StructuralElement
 * @see ThermalElement
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../core/config.h"
#include "../core/parallel.h"
#include "element/element_structural.h"
#include "element/element_thermal.h"
#include "model.h"

#include <utility>
#include <vector>

namespace fem::model {

/**
 * Recovers unaveraged element-nodal stress and strain.
 *
 * Every supporting structural element supplies its natural nodal recovery
 * coordinates and writes six-component stress and strain values into
 * `ELEMENT_NODAL` storage. The compiled element-nodal offsets keep parallel
 * writes disjoint. Unwritten element-local rows remain NaN and are omitted
 * by Field::to_nodal() when a global nodal result is requested.
 *
 * Elements without nodal recovery coordinates do not participate.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @param linearization Displacement expansion state; null denotes zero.
 * @return Element-nodal stress and strain fields.
 */
std::tuple<Field, Field> Model::compute_stress(
    Field&       target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const Index element_count = static_cast<Index>(_data->elements.size());

    // Allocate disjoint element-nodal recovery fields
    Field element_stress = _data->create_field_("ELEMENT_NODAL_STRESS", FieldDomain::ELEMENT_NODAL, 6);
    Field element_strain = _data->create_field_("ELEMENT_NODAL_STRAIN", FieldDomain::ELEMENT_NODAL, 6);

    // Recover element-nodal values into non-overlapping compiled row ranges
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!element) return;

            if (auto sel = element->as<StructuralElement>()) {
                RowMatrix rst = sel->stress_strain_nodal_rst();
                if (rst.rows() == 0) return;

                logging::error(rst.rows() == sel->n_nodes(),
                    "Element ", sel->elem_id, " returned ", rst.rows(),
                    " nodal stress coordinates, expected ", sel->n_nodes());
                sel->compute_stress_strain(
                    &element_strain,
                    &element_stress,
                    target_displacement,
                    target_temperature,
                    rst,
                    base_displacement,
                    base_temperature
                );
            }
        }, Index(8));

    return {std::move(element_stress), std::move(element_strain)};
}

/**
 * Recovers accumulated equivalent plastic strain at element nodes.
 *
 * Supporting structural elements read PEEQ directly from accepted constitutive
 * history and write one scalar value per element node. Solid and shell
 * formulations perform their integration-point extrapolation before this
 * model-level operation. The compiled element-nodal offsets keep those writes
 * disjoint until a nodal projection is explicitly requested.
 *
 * Elements without a PEEQ-capable material or formulation retain NaN and do
 * not dilute the nodal average.
 *
 * @return Scalar ELEMENT_NODAL field named PEEQ.
 */
Field Model::compute_peeq_element_nodal() {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const Index element_count = static_cast<Index>(_data->elements.size());

    // Unsupported element-local PEEQ rows remain unwritten (NaN)
    Field element_peeq = _data->create_field_("ELEMENT_NODAL_PEEQ", FieldDomain::ELEMENT_NODAL, 1);

    // Let each structural formulation recover accepted material history locally
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!element) return;

            if (auto sel = element->as<StructuralElement>()) {
                sel->compute_peeq(element_peeq);
            }
        }, Index(8));

    // Keep unaveraged material-history results for optional nodal projection
    return element_peeq;
}

/**
 * Retains the existing direct nodal PEEQ API for model callers.
 */
Field Model::compute_peeq_nodal() {
    return compute_peeq_element_nodal().to_nodal(*_data, "PEEQ");
}

/**
 * Recovers element-nodal stress on the top and bottom shell faces.
 *
 * Shell recovery coordinates use the formulation-provided natural nodal points
 * with the through-thickness coordinate set to `+1` for the top face and `-1`
 * for the bottom face. Non-shell structural elements use their unchanged nodal
 * recovery coordinates for both fields. Unwritten rows remain NaN and are
 * omitted by Field::to_nodal() when nodal results are requested.
 *
 * @param displacement Global nodal displacement field used for recovery.
 * @param linearization Displacement expansion state; null denotes zero.
 * @return Element-nodal top and bottom stress fields.
 */
std::tuple<Field, Field> Model::compute_stress_top_bot(
    Field&       target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const Index element_count = static_cast<Index>(_data->elements.size());

    // Allocate separate element-nodal fields for both faces
    Field element_top = _data->create_field_("ELEMENT_NODAL_STRESS_TOP", FieldDomain::ELEMENT_NODAL, 6);
    Field element_bot = _data->create_field_("ELEMENT_NODAL_STRESS_BOT", FieldDomain::ELEMENT_NODAL, 6);

    // Evaluate both faces in each element's disjoint compiled row range
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!element) return;

            if (auto sel = element->as<StructuralElement>()) {
                const RowMatrix base_rst = sel->stress_strain_nodal_rst();
                if (base_rst.rows() == 0) return;

                RowMatrix rst_bot = base_rst;
                RowMatrix rst_top = base_rst;

                // Select the two natural thickness faces for shell formulations
                if (sel->is_shell()) {
                    for (int i = 0; i < base_rst.rows(); ++i) {
                        rst_bot(i, 2) = -1;
                        rst_top(i, 2) =  1;
                    }
                }

                sel->compute_stress_strain(
                    nullptr, &element_bot,
                    target_displacement, target_temperature,
                    rst_bot,
                    base_displacement, base_temperature);
                sel->compute_stress_strain(
                    nullptr, &element_top,
                    target_displacement, target_temperature,
                    rst_top,
                    base_displacement, base_temperature);
            }
        }, Index(8));

    return {std::move(element_top), std::move(element_bot)};
}

/**
 * Recovers shell forces and moments in the unaveraged element-nodal domain.
 *
 * Each shell writes its eight generalized resultants into its own compiled
 * row range, so adjacent elements may preserve discontinuities and be evaluated
 * concurrently. A NODE output request projects these values separately.
 *
 * @param target_displacement Global nodal displacement field.
 * @param target_temperature Target temperature field.
 * @param base_displacement Linearization displacement field.
 * @param base_temperature Linearization temperature field.
 * @return Eight-component ELEMENT_NODAL result field.
 */
Field Model::compute_shell_resultants(
    Field&       target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    logging::error(_data->element_nodal_offsets != nullptr,
        "Model: element-nodal offsets are required for shell resultants");

    const Index element_count = static_cast<Index>(_data->elements.size());
    Field resultants = _data->create_field_(
        "ELEMENT_NODAL_SHELL_RESULTANTS", FieldDomain::ELEMENT_NODAL, 8);

    // Each element writes a disjoint segment; no global nodal reduction is needed
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!element) return;
            if (auto sel = element->as<StructuralElement>()) {
                sel->compute_shell_section_forces(
                    resultants, target_displacement, target_temperature,
                    base_displacement, base_temperature);
            }
        }, Index(8));

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
    Field compliance = _data->create_field_("COMPLIANCE", FieldDomain::ELEMENT, 1, false);
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
    Field results = _data->create_field_("ORIENTATION_GRAD", FieldDomain::ELEMENT, 3, false);
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
    Field volumes = _data->create_field_("VOLUME", FieldDomain::ELEMENT, 1, false);
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
 * @param target_displacement Requested nodal displacement state.
 * @param target_temperature Requested nodal temperature state.
 * @param base_displacement Base displacement state; null denotes zero.
 * @param base_temperature Base temperature state; null denotes stress-free temperature.
 * @return Six-component element-nodal beam section-force field.
 */
Field Model::compute_section_forces(
    Field&       target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    // Validate and access the compiled element-nodal enumeration
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate zeroed storage for every compiled element-node row
    Field beam_forces = _data->create_field_("BEAM_SECTION_FORCES", FieldDomain::ELEMENT_NODAL, 6, false);
    beam_forces.set_zero();

    // Recover section forces into each element's disjoint nodal row range
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!el) return;
            if (auto sel = el->as<StructuralElement>()) {
                sel->compute_beam_section_forces(
                    beam_forces,
                    target_displacement,
                    target_temperature,
                    base_displacement,
                    base_temperature);
            }
        }, Index(32));

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

    const Index element_count       = static_cast<Index>(_data->elements.size());

    // Allocate zeroed storage for every compiled element-node row
    Field shear_flow = _data->create_field_("SHEAR_FLOW", FieldDomain::ELEMENT_NODAL, 1, false);
    shear_flow.set_zero();

    // Recover shear flow into each element's disjoint nodal row range
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!el) return;
            if (auto sel = el->as<StructuralElement>()) {
                sel->compute_shear_flow(shear_flow, displacement);
            }
        }, Index(32));

    return shear_flow;
}

/**
 * Recovers conductive heat flux as an element-nodal result field.
 *
 * Each thermal element first produces one heat-flux vector per local node in its
 * disjoint `ELEMENT_NODAL` row range. How those values are obtained is owned by
 * the concrete formulation: solid elements evaluate Fourier's law at their
 * integration points and apply their topology-specific extrapolation matrix,
 * while other thermal formulations may evaluate the gradient directly at their
 * nodes when that is their natural recovery scheme.
 *
 * Element recovery is parallel because no two elements write the same
 * `ELEMENT_NODAL` rows. As with stress recovery, the caller decides whether
 * to preserve these discontinuous values or project them onto shared nodes
 * using `Field::to_nodal()`.
 *
 * No model-wide integration-point heat-flux field is created or exported.
 *
 * @param temperature Scalar global nodal temperature field.
 * @return Three-component `ELEMENT_NODAL` heat-flux field.
 */
Field Model::compute_heat_flux(const Field& temperature) {
    // Validate the scalar primary field and compiled element-nodal enumeration
    // required by the parallel local recovery path.
    logging::error(temperature.domain == FieldDomain::NODE,
        "Model: temperature field must use the NODE domain");
    logging::error(temperature.components == 1,
        "Model: temperature field must have exactly one component");
    logging::error(temperature.rows == _data->field_rows(FieldDomain::NODE),
        "Model: temperature field does not match the compiled nodal domain");
    logging::error(_data->element_nodal_offsets != nullptr,
        "Model: element nodal offsets are not initialized");

    const Index element_count = static_cast<Index>(_data->elements.size());

    // Thermal elements write their values; unsupported rows remain NaN
    Field element_heat_flux = _data->create_field_("ELEMENT_NODAL_HEAT_FLUX", FieldDomain::ELEMENT_NODAL, 3);

    // Recover each thermal element into its disjoint element-nodal field range
    parallel::for_index(element_count, global_config.max_threads,
        [&](Index elem_idx, int /*worker*/) {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (!element) {
                return;
            }

            if (auto* thermal = element->as<ThermalElement>()) {
                thermal->compute_heat_flux(element_heat_flux, temperature);
            }
        }, Index(8));


    // Preserve the discontinuous element-local values for output selection
    return element_heat_flux;
}


}
