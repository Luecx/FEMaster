//
// Created by f_eggers on 10.02.2026.
//

#ifndef COMPILE_SH_INERTIA_RELIEF_H
#define COMPILE_SH_INERTIA_RELIEF_H

#include "../../model/model_data.h"

namespace fem {

void apply_inertia_relief(model::ModelData& model_data,
                          model::Field& global_load_mat,
                          bool consider_point_masses = true);
}

#endif    // COMPILE_SH_INERTIA_RELIEF_H
