//
// Created by f_eggers on 10.02.2026.
//

#ifndef COMPILE_SH_INERTIA_RELIEF_H
#define COMPILE_SH_INERTIA_RELIEF_H

#include "../../model/model_data.h"

namespace fem {

void apply_inertia_relief(model::ModelData& model_data, model::Field& global_load_mat);

}

#endif    // COMPILE_SH_INERTIA_RELIEF_H
