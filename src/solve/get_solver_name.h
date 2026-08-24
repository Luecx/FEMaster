#pragma once

#include "solve_device.h"
#include "solve_method.h"
#include "sparse/solve_sparse_direct.h"

#include <string>

namespace fem::solver {

std::string get_solver_name(SolverDevice device,
                            SolverMethod method,
                            DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD);

} // namespace fem::solver
