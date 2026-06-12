/**
 * @file nonlinear_static.h
 * @brief Declares the nonlinear static load case.
 */

#pragma once

#include "loadcase.h"
#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>
#include <vector>

namespace fem {
namespace loadcase {

/**
 * @struct NonlinearStatic
 * @brief Incremental quasi-Newton updated-Lagrangian static analysis.
 */
struct NonlinearStatic : public LoadCase {
    std::vector<std::string> supps;
    std::vector<std::string> loads;
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;
    std::string stiffness_file;

    int num_increments = 10;
    int max_iterations = 20;
    Precision tolerance = Precision(1e-8);
    bool regularize_zero_stiffness_rows = true;
    Precision zero_stiffness_regularization_alpha = Precision(1e-4);

    NonlinearStatic(ID id, reader::Writer* writer, model::Model* model);

    void run() override;
};

} // namespace loadcase
} // namespace fem
