#include "beam_generalized_strain.h"

namespace fem {

BeamGeneralizedStrain::BeamGeneralizedStrain(const Vec6& values)
    : values_(values) {}

const Vec6& BeamGeneralizedStrain::values() const {
    return values_;
}

Vec6& BeamGeneralizedStrain::values() {
    return values_;
}

} // namespace fem
