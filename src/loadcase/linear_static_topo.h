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
    // mesh independence filtering
    Precision filter_radius = 0;
    Precision gaussian_sigma = 0;

    private:
    ElementData filter(ElementData& data);

    public:
    void run() override;
};

}
}