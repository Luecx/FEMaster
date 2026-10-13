#pragma once

#include "../../core/types_eig.h"
#include "../../model/model_data.h"
#include "../../section/pretension_section.h"
#include "equation.h"

namespace fem::constraint {

Equations get_pretension_equations(
    SystemDofIds& system_dof_ids,
    model::ModelData& model_data,
    pretension::PretensionSection& section);

}
