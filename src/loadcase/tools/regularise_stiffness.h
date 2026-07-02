#pragma once

#include "../../core/types_eig.h"

namespace fem {

int regularise_stiffness(SparseMatrix& matrix, Precision alpha);

} // namespace fem
