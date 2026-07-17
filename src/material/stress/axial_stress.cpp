#include "axial_stress.h"

namespace fem {

AxialStress::AxialStress(Precision value)
    : value_(value) {}

Precision AxialStress::value() const {
    return value_;
}

Precision& AxialStress::value() {
    return value_;
}

} // namespace fem
