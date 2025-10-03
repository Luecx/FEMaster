#include "loadcase.h"

#include "../constraints/constraint_report.h"
#include "../core/logging.h"

namespace fem {
namespace loadcase {

void LoadCase::report_constraint_groups(const model::Model::ConstraintGroups& groups) const {
    using fem::logging::info;
    if (!report_constraints) {
        return;
    }

    const constraint::EquationFormatOptions fmt{};

    const struct {
        const char* name;
        const constraint::Equations& eqs;
    } categories[] = {
        {"Supports",   groups.supports},
        {"Connectors", groups.connectors},
        {"Couplings",  groups.couplings},
        {"Ties",       groups.ties},
        {"Other",      groups.others},
    };

    info(true, "Constraint overview for loadcase ", m_id, " (", m_model->element_dims, "D elements)");
    logging::up();
    for (const auto& category : categories) {
        info(true, category.name, " (", category.eqs.size(), ")");
        if (!category.eqs.empty()) {
            logging::up();
            for (const auto& line : constraint::format_equations(category.eqs, fmt)) {
                info(true, "    ", line);
            }
            logging::down();
        }
    }
    logging::down();
}

} // namespace loadcase
} // namespace fem
