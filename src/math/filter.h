#pragma once

#include "../core/types.h"
#include "../solve/solver.h"

namespace fem {
namespace math {

DynamicVector filter_gaussian(DynamicMatrix &coords, DynamicVector &values, Precision radius, Precision sigma, solver::SolverDevice device);

}
}    // namespace fem