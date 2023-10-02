#pragma once

#include "loadcase.h"
#include "../solve/solver.h"
#include "../math/interpolate.h"

namespace fem{
namespace loadcase{

#define TIMED_EXECUTION(function_call, description) \
    do { \
        timer.start(); \
        function_call; \
        timer.stop(); \
        logging::info(true, "It took", timer.elapsed(), "milliseconds to run", description); \
    } while(0)

struct LinearStatic : public LoadCase{

    explicit LinearStatic(ID id, reader::Writer* writer, model::Model* model);

    std::vector<std::string> supps;
    std::vector<std::string> loads;

    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::INDIRECT;

    public:
    virtual void run() override;


};

#undef TIMED_EXECUTION

}
}