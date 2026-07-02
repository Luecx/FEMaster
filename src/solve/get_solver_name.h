#pragma once

#include "solve_device.h"
#include "solve_method.h"

#include <string>

namespace fem::solver {

std::string get_solver_name(SolverDevice device, SolverMethod method);

} // namespace fem::solver
