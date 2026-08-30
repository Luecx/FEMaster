/**
 * @file model_field_mapping.cpp
 * @brief Implements model-level element weights for element-nodal projection.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "element/element_structural.h"
#include "element/element_thermal.h"
#include "model.h"

namespace fem::model {

/**
 * Builds unit element weights for element-nodal to nodal field projection.
 *
 * Structural and thermal capabilities can be selected independently. An element
 * receives weight one when it implements at least one requested capability and
 * zero otherwise. Elements implementing both capabilities still receive a
 * single unit weight.
 *
 * @param structural Include elements implementing `StructuralElement`.
 * @param thermal Include elements implementing `ThermalElement`.
 * @return Scalar `ELEMENT` field containing zero or one for every dense element.
 */
Field Model::build_field_mapping_weights(bool structural, bool thermal) const {
    logging::error(structural || thermal,
        "Model: field mapping weights require at least one element capability");

    const Index element_count = static_cast<Index>(_data->elements.size());
    Field weights{"FIELD_MAPPING_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    weights.set_zero();

    for (const auto& element : _data->elements) {
        if (element == nullptr) continue;

        const bool selected =
            (structural && element->as<StructuralElement>() != nullptr)
         || (thermal    && element->as<ThermalElement>()    != nullptr);
        if (!selected) continue;

        logging::error(element->elem_id >= 0
                    && static_cast<Index>(element->elem_id) < element_count,
            "Model: element id out of range while building field mapping weights: ",
            element->elem_id);
        weights(static_cast<Index>(element->elem_id), 0) = Precision(1);
    }

    return weights;
}

} // namespace fem::model
