#include "axial_strain.h"

namespace fem {

AxialStrain::AxialStrain(Precision value)
    : value_(value) {}

Precision AxialStrain::value() const {
    return value_;
}

Precision& AxialStrain::value() {
    return value_;
}

} // namespace fem
