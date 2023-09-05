#pragma once

#include "loadcase.h"
#include "linear_static.h"
#include "../solve/solver.h"

namespace fem{
namespace loadcase{

#define TIMED_EXECUTION(function_call, description) \
    do { \
        timer.start(); \
        function_call; \
        timer.stop(); \
        logging::info(true, "It took", timer.elapsed(), "milliseconds to run", description); \
    } while(0)

struct LinearStaticTopo : public LinearStatic{

    explicit LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model);

    ElementData density;
    Precision exponent = 1;

    public:
    void run() override;
};

#undef TIMED_EXECUTION

}
}