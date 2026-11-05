#include "pretension.h"

#include "../../constraints/types/pretension.h"
#include "../../core/logging.h"
#include "../../data/field.h"

#include <utility>

namespace fem::bc {

PretensionCondition::PretensionCondition(pretension::PretensionSection::Ptr section_,
                                         Mode mode_, Precision force)
    : section(std::move(section_)), mode(mode_),
      value(mode == Mode::Hold && section ? section->last_solved_gap : force) {
    logging::error(section != nullptr, "PRETENSION: missing section");
    logging::error(mode != Mode::Hold || section->has_solved_gap,
        "PRETENSION HOLD: section ", section->name,
        " has no successfully completed prior step");
    logging::error(std::isfinite(value), "PRETENSION: prescribed value must be finite");
}

void PretensionCondition::apply(model::ModelData&, model::Field& rhs,
                                constraint::Equations& equations,
                                const SystemDofIds& dofs, TripletList&, Precision scale) {
    // Constraint collection supplies numbering and no RHS; load assembly
    // supplies an RHS and no numbering, as for the other physical conditions.
    if (dofs.rows() != 0) {
        auto rows = constraint::get_pretension_equations(
            dofs, *section, mode == Mode::Force, mode == Mode::Hold ? value : Precision(0));
        equations.insert(equations.end(), rows.begin(), rows.end());
    }
    if (rhs.rows == 0 || mode != Mode::Force) return;
    logging::error(!section->interface_pairs.empty(), "PRETENSION: section is not prepared");
    const Vec3 force = scale * value / static_cast<Precision>(section->interface_pairs.size())
        * section->axis_direction;
    for (const auto& pair : section->interface_pairs)
        for (Dim c = 0; c < 3; ++c) {
            rhs(pair.side_a, c) -= force(c);
            rhs(pair.side_b, c) += force(c);
        }
}

std::string PretensionCondition::str() const {
    return "PRETENSION " + section->name + (mode == Mode::Force ? " FORCE " : " HOLD ")
        + std::to_string(value);
}

}
