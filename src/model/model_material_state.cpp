/**
 * @file model_material_state.cpp
 * @brief Implements model-wide material-point state initialization.
 *
 * Constitutive history is stored in one `ELEMENT_MP` field whose rows follow
 * the material-point enumeration created during model compilation. Every row
 * reserves a uniform component width equal to the largest state vector required
 * by any assigned elastic law in the compiled model.
 *
 * This file only defines the initial value of that storage. Trial, committed and
 * rollback semantics during nonlinear solution remain the responsibility of the
 * nonlinear state manager and the constitutive evaluation workflow.
 *
 * @see Model
 * @see material::Elasticity
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "model.h"

#include "../core/logging.h"
#include "../material/elasticity.h"

#include <algorithm>

namespace fem::model {

/**
 * Determines the uniform component width required by material history fields.
 *
 * The routine scans all compiled elements with an assigned section, material
 * and elastic law. Stateless laws contribute zero components; elements without
 * a complete constitutive assignment are ignored. The maximum is suitable for
 * every row of the model-wide `ELEMENT_MP` state field.
 *
 * @return Largest constitutive state-vector size used by any compiled element.
 */
Index Model::maximum_material_state_size() const {
    Index maximum_size = 0;

    // Accumulate the largest state requirement across assigned elastic laws
    for (const auto& element : _data->elements) {
        if (!element || !element->_section || !element->_section->material_) {
            continue;
        }

        const auto& material = element->_section->material_;
        if (!material->has_elasticity()) {
            continue;
        }

        maximum_size = std::max(maximum_size, material->elasticity()->state_size());
    }

    return maximum_size;
}

/**
 * Initializes every enumerated material-point state in a compatible field.
 *
 * The field must use the `ELEMENT_MP` domain, contain exactly the row count
 * established by compiled material-point enumeration and provide at least the
 * maximum component width required by the model. All storage is cleared first.
 * Each stateful elastic law then initializes the contiguous state vector at
 * every integration-point/material-point row owned by its element.
 *
 * Elements without a section, material, elastic law or state variables leave
 * their allocated rows at zero.
 *
 * @param state Model-wide material-point history field to initialize in place.
 */
void Model::initialize_material_state(Field& state) const {
    // Validate the field domain and compiled material-point dimensions
    logging::error(state.domain == FieldDomain::ELEMENT_MP,
        "Material state field must use ELEMENT_MP domain");

    const Index material_points = _data->field_rows(FieldDomain::ELEMENT_MP);
    logging::error(state.rows == material_points,
        "Material state field has ", state.rows,
        " rows, expected ", material_points);

    const Index maximum_size = maximum_material_state_size();
    logging::error(state.components >= maximum_size,
        "Material state field has ", state.components,
        " components, expected at least ", maximum_size);

    // Establish a deterministic zero state for unused rows and components
    state.set_zero();

    // Initialize the constitutive state vector at every represented material point
    for (const auto& element : _data->elements) {
        if (!element || !element->_section || !element->_section->material_) {
            continue;
        }

        const auto& material = element->_section->material_;
        if (!material->has_elasticity()) {
            continue;
        }

        const auto& elasticity = material->elasticity();
        const Index state_size = elasticity->state_size();
        if (state_size == 0) {
            continue;
        }

        logging::error(state_size <= state.components,
            "Material state size ", state_size,
            " exceeds field width ", state.components,
            " for element ", element->elem_id);

        // Resolve each element-local integration/material point to its global row
        for (Index ip = 0; ip < static_cast<Index>(element->num_ip()); ++ip) {
            for (Index mp = 0; mp < element->num_mp_per_ip(); ++mp) {
                const Index row = element->mp_index(ip, mp);
                logging::error(row < state.rows,
                    "Material point row ", row,
                    " is out of range for element ", element->elem_id);

                elasticity->initialize_state(state.data() + row * state.components);
            }
        }
    }
}

} // namespace fem::model
