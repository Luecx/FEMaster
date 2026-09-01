/**
 * @file model_field_mapping.cpp
 * @brief Implements model-level element weights for element-nodal projection.
 *
 * Element-nodal result recovery first computes quantities independently in every
 * element's disjoint output range and subsequently projects those values to the
 * shared global nodes. The scalar weights built here select which element
 * capabilities participate in that averaging step.
 *
 * Structural and thermal recovery may use the same compiled topology but must not
 * average through elements that do not provide the requested physical quantity.
 * The capability flags therefore allow stress recovery to select structural
 * elements only and heat-flux recovery to select thermal elements only.
 *
 * @see Model::build_field_mapping_weights
 * @see ModelData::element_nodal_to_nodal
 * @see StructuralElement
 * @see ThermalElement
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "element/element_structural.h"
#include "element/element_thermal.h"
#include "model.h"

namespace fem::model {

/**
 * Builds unit participation weights for element-nodal to nodal field projection.
 *
 * Every dense compiled element receives one scalar `ELEMENT` weight. The weight
 * is one when the element implements at least one requested capability and zero
 * otherwise. An element implementing both structural and thermal behavior still
 * receives one rather than two because the field is a participation mask, not a
 * multiplicity measure.
 *
 * These weights are consumed by `ModelData::element_nodal_to_nodal()` to exclude
 * unrelated elements from nodal averaging. For example, thermal heat flux should
 * not be diluted by a structural-only element sharing the same node.
 *
 * @param structural Include elements implementing `StructuralElement`.
 * @param thermal Include elements implementing `ThermalElement`.
 * @return Scalar `ELEMENT` field containing zero or one for every dense element.
 */
Field Model::build_field_mapping_weights(bool structural, bool thermal) const {
    // At least one physical capability must be selected; an all-zero mapping has
    // no meaningful recovery interpretation and likely indicates a caller error.
    logging::error(structural || thermal,
        "Model: field mapping weights require at least one element capability");

    // Allocate one scalar participation value per dense compiled element. Zero is
    // the correct default for empty slots and elements outside the requested
    // physical capability set.
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field weights{"FIELD_MAPPING_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    weights.set_zero();

    // Classify each initialized element by its implemented physical capabilities.
    for (const auto& element : _data->elements) {
        if (element == nullptr) continue;

        // Select the element when any requested capability is implemented. The
        // boolean combination deliberately collapses dual-capability elements to
        // one unit participation weight.
        const bool selected =
            (structural && element->as<StructuralElement>() != nullptr)
         || (thermal    && element->as<ThermalElement>()    != nullptr);
        if (!selected) continue;

        // Compiled element identifiers index the dense ELEMENT field directly and
        // must therefore lie within the allocated mapping range.
        logging::error(element->elem_id >= 0
                    && static_cast<Index>(element->elem_id) < element_count,
            "Model: element id out of range while building field mapping weights: ",
            element->elem_id);

        weights(static_cast<Index>(element->elem_id), 0) = Precision(1);
    }

    return weights;
}

} // namespace fem::model
