/**
 * @file model_material_state.cpp
 * @brief Implements model-wide material-point state initialization.
 *
 * Material history variables are stored in one dense `ELEMENT_MP` field. Every
 * row belongs to one globally enumerated material point and every row has the
 * maximum state width required by any constitutive law that is actually assigned
 * to an element. Materials use only the leading
 * `ConstitutiveLaw::state_size()` entries of their rows.
 *
 * The model initializes storage but does not own nonlinear trial buffers. The
 * active nonlinear load case owns committed and trial fields and publishes them
 * through `ModelData::material_state_old` and `material_state_new`.
 *
 * @see Model::maximum_material_state_size
 * @see Model::initialize_material_state
 * @see material::ConstitutiveLaw::state_size
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "model.h"

#include "../core/logging.h"

#include <algorithm>

namespace fem {
namespace model {

/**
 * Determines the common component width required by nonlinear material-state
 * fields.
 *
 * Only materials assigned through an element section participate. Elements
 * without a section, without a material or without elasticity are skipped. A
 * single dense field uses the maximum state width so every globally enumerated
 * material point has a constant row stride; individual materials consume only
 * their leading `ConstitutiveLaw::state_size()` components.
 *
 * @return Maximum constitutive state size of all assigned materials, or zero
 *         when the model is entirely stateless.
 */
Index Model::maximum_material_state_size() const {
    Index maximum_size = 0;

    for (const auto& element : _data->elements) {
        if (!element || !element->_section || !element->_section->material_) {
            continue;
        }

        const auto& material = element->_section->material_;
        if (!material->has_elasticity()) {
            continue;
        }

        maximum_size = std::max(
            maximum_size,
            material->constitutive_law().state_size()
        );
    }

    return maximum_size;
}

/**
 * Initializes every enumerated material-point row to its constitutive reference
 * history.
 *
 * The supplied field must cover the complete `ELEMENT_MP` domain and provide at
 * least the maximum state width required by any assigned constitutive law. The
 * field is zeroed first so unused trailing components and stateless rows have a
 * stable value. Stateful laws then initialize each integration-point/material-
 * point row belonging to their element.
 *
 * Element offsets must already have been established by
 * `ModelData::initialize_element_enumeration()`. The operation initializes
 * storage only; it neither binds the field as active state nor commits nonlinear
 * history.
 *
 * @param state Material-state field to initialize in global `ELEMENT_MP` order.
 */
void Model::initialize_material_state(Field& state) const {
    logging::error(state.domain == FieldDomain::ELEMENT_MP,
        "Material state field must use ELEMENT_MP domain");
    logging::error(state.rows == static_cast<Index>(_data->max_material_points),
        "Material state field has ", state.rows,
        " rows, expected ", _data->max_material_points);

    const Index maximum_size = maximum_material_state_size();
    logging::error(state.components >= maximum_size,
        "Material state field has ", state.components,
        " components, expected at least ", maximum_size);

    state.set_zero();

    for (const auto& element : _data->elements) {
        if (!element || !element->_section || !element->_section->material_) {
            continue;
        }

        const auto& material = element->_section->material_;
        if (!material->has_elasticity()) {
            continue;
        }

        const auto& law = material->constitutive_law();
        const Index state_size = law.state_size();
        if (state_size == 0) {
            continue;
        }

        logging::error(state_size <= state.components,
            "Material state size ", state_size,
            " exceeds field width ", state.components,
            " for element ", element->elem_id);

        for (Index ip = 0; ip < static_cast<Index>(element->num_ip()); ++ip) {
            for (Index mp = 0; mp < element->num_mp_per_ip(); ++mp) {
                const Index row = element->mp_index(ip, mp);
                logging::error(row >= 0 && row < state.rows,
                    "Material point row ", row,
                    " is out of range for element ", element->elem_id);

                law.initialize_state(
                    state.data() + row * state.components
                );
            }
        }
    }
}

} // namespace model
} // namespace fem
