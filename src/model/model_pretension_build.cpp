#include "model.h"

namespace fem::model {

Field Model::build_pretension_force_matrix() {
    Field force{"PRETENSION_FORCE", FieldDomain::NODE,
                _data->field_rows(FieldDomain::NODE), 3};
    force.set_zero();
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->nodal_control || section->state != pretension::State::Loading ||
            section->control != pretension::Control::Force ||
            section->interface_pairs.empty()) continue;
        force(section->control_node, 0) = section->prescribed_value;
    }
    return force;
}

Field Model::build_pretension_gap_matrix(const Field& displacement) {
    logging::error(displacement.domain == FieldDomain::NODE && displacement.components >= 3,
        "Pretension gap output requires a nodal displacement field");
    Field gap{"PRETENSION_GAP", FieldDomain::NODE,
              _data->field_rows(FieldDomain::NODE), 1};
    gap.set_zero();
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->interface_pairs.empty()) continue;
        const Vec3 axis = section->axis_direction.normalized();
        for (const auto& pair : section->interface_pairs) {
            const Vec3 relative =
                displacement.row_vec3(static_cast<Index>(pair.side_b)) -
                displacement.row_vec3(static_cast<Index>(pair.side_a));
            const Precision measured = axis.dot(relative);
            gap(static_cast<Index>(pair.side_a), 0) = measured;
            gap(static_cast<Index>(pair.side_b), 0) = measured;
        }
    }
    return gap;
}

void Model::capture_pretension_gaps(const Field& displacement) {
    logging::error(displacement.domain == FieldDomain::NODE && displacement.components >= 3,
        "Capturing pretension gaps requires a nodal displacement field");
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->state == pretension::State::Open ||
            section->interface_pairs.empty()) continue;
        section->last_solved_gap = displacement(section->control_node, 0);
        section->has_solved_gap = true;
        logging::info(true, "PRETENSION '", section->name,
            "': solved mean gap = ", section->last_solved_gap);
    }
}

} // namespace fem::model
