/**
 * @file model_material_state.cpp
 * @brief Implements model-wide material-point state initialization.
 *
 * Material history variables are stored in one dense `ELEMENT_MP` field. Every
 * row belongs to one globally enumerated material point and every row has the
 * maximum state width required by any material that is actually assigned to an
 * element. Materials use only the leading `Elasticity::state_size()` entries of
 * their rows.
 *
 * The model initializes storage but does not own nonlinear trial buffers. The
 * active load case owns committed and trial fields and binds only the currently
 * active state through `ModelData::material_state`.
 *
 * @see Model::maximum_material_state_size
 * @see Model::initialize_material_state
 * @see material::Elasticity::state_size
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "model.h"

#include "../core/logging.h"
#include "../material/elasticity.h"

#include <algorithm>

namespace fem {
namespace model {

/**
 * Determines the common component width required by the active material-state
 * field.
 *
 * Only materials assigned through an element section participate. Elements
 * without a section, without a material or without elasticity are skipped. A
 * single dense field uses the maximum state width so every globally enumerated
 * material point has a constant row stride; individual materials consume only
 * their leading `Elasticity::state_size()` components.
 *
 * @return Maximum constitutive state size of all assigned elasticities, or zero
 *         when the model is entirely stateless.
 */
Index Model::maximum_material_state_size() const {
    Index maximum_size = 0;

    // Inspect only constitutive models reachable from active element sections
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
            material->elasticity()->state_size()
        );
    }

    return maximum_size;
}

/**
 * Initializes every enumerated material-point row to its constitutive reference
 * history.
 *
 * The supplied field must cover the complete `ELEMENT_MP` domain and provide at
 * least the maximum state width required by any assigned elasticity. The field
 * is zeroed first so unused trailing components and stateless rows have a stable
 * value. Stateful materials then initialize each integration-point/material-
 * point row belonging to their element through the common elasticity API.
 *
 * Element offsets must already have been established by
 * `ModelData::initialize_element_enumeration()`. The operation initializes
 * storage only; it neither binds the field as active state nor commits nonlinear
 * history.
 *
 * @param state Material-state field to initialize in global `ELEMENT_MP` order.
 */
void Model::initialize_material_state(Field& state) const {
    // Validate the global material-point layout before touching the field
    logging::error(state.domain == FieldDomain::ELEMENT_MP,
        "Material state field must use ELEMENT_MP domain");
    logging::error(state.rows == static_cast<Index>(_data->max_material_points),
        "Material state field has ", state.rows,
        " rows, expected ", _data->max_material_points);

    const Index maximum_size = maximum_material_state_size();
    logging::error(state.components >= maximum_size,
        "Material state field has ", state.components,
        " components, expected at least ", maximum_size);

    // Establish deterministic values for stateless points and unused trailing
    // components before material-specific initialization
    state.set_zero();

    // Initialize every stateful material point in element, integration-point
    // and section-material-point order
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

        // The element owns the mapping from its local IP/MP pair into the
        // globally flattened material-point field
        for (Index ip = 0; ip < static_cast<Index>(element->num_ip()); ++ip) {
            for (Index mp = 0; mp < element->num_mp_per_ip(); ++mp) {
                const Index row = element->mp_index(ip, mp);
                logging::error(row >= 0 && row < state.rows,
                    "Material point row ", row,
                    " is out of range for element ", element->elem_id);

                elasticity->initialize_state(
                    state.data() + row * state.components
                );
            }
        }
    }
}

} // namespace model
} // namespace fem
