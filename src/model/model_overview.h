//
// Model overview pretty-printer
//

#pragma once

#include "model.h"

namespace fem { namespace model {

// Emits a structured overview of the model to the logger
void print_model_overview(const Model& mdl);

} } // namespace fem::model

