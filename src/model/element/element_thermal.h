
#pragma once

#include "element.h"
#include <functional>

namespace fem {
namespace model {

struct ThermalElement {
    virtual ~ThermalElement() = default;

    // conductivity and capacity matrix for head transfer problems
    virtual MapMatrix conductivity(Precision* buffer) = 0;
    virtual MapMatrix capacity    (Precision* buffer) = 0;

    virtual void compute_heat_flux(
        Field&       heat_flux,
        const Field& temperature
    ) = 0;
};


} // namespace model
} // namespace fem