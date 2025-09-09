// loadcases/linear_buckling.h
#pragma once
#include "loadcase.h"
#include "../solve/solver.h"

namespace fem { namespace loadcase {

struct LinearBuckling : public LoadCase {
    explicit LinearBuckling(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues);

    std::vector<std::string> supps;
    std::vector<std::string> loads;
    int num_eigenvalues;

    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;

    virtual void run() override;
};

}} // ns
