#pragma once

#include "../condition.h"
#include "../../section/pretension_section.h"

namespace fem::bc {

// Immutable step prescription. Accepted solution history stays on the section.
struct PretensionCondition : Condition {
    enum class Mode { Force, Hold };

    const pretension::PretensionSection::Ptr section;
    const Mode mode;
    const Precision value; // Total FORCE or the accepted gap frozen for HOLD.

    PretensionCondition(pretension::PretensionSection::Ptr section,
                        Mode mode, Precision force = Precision(0));

    void apply(model::ModelData& data, model::Field& rhs,
               constraint::Equations& equations, const SystemDofIds& dofs,
               TripletList& matrix, Precision scale = Precision(1)) override;
    std::string str() const override;
};

}
