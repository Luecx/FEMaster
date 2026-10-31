#pragma once

#include "../../core/types_eig.h"
#include "../../section/pretension_section.h"
#include "equation.h"

namespace fem::constraint {

Equations get_pretension_equations(
    const SystemDofIds& system_dof_ids,
    const pretension::PretensionSection& section,
    bool force_control, Precision normal_gap);

}
