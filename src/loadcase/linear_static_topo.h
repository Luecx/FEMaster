#pragma once

#include "loadcase.h"
#include "linear_static.h"
#include "../solve/solver.h"

namespace fem{
namespace loadcase{

struct LinearStaticTopo : public LinearStatic{

    explicit LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model);

    ElementData density;
    Precision exponent = 1;

    public:
    void run() override;
};

}
}