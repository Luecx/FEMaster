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
 * active load case creates and binds committed and trial fields through the
 * semantic pointers in `ModelData`.
 *
 * @see Model::maximum_material_state_size
 * @see Model::initialize_material_state
 * @see material::Elasticity::state_size
 */

#include "model.h"

#include "../core/logging.h"
#include "../material/elasticity.h"

#include <algorithm>

namespace fem {
namespace model {

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
            material->elasticity()->state_size()
        );
    }

    return maximum_size;
}

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

        const auto& elasticity = material->elasticity();
        const Index state_size = elasticity->state_size();
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

                elasticity->initialize_state(
                    state.data() + row * state.components
                );
            }
        }
    }
}

} // namespace model
} // namespace fem
