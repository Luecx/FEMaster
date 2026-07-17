#include "beam_stress_resultants.h"

namespace fem {

BeamStressResultants::BeamStressResultants(const Vec6& values)
    : values_(values) {}

const Vec6& BeamStressResultants::values() const {
    return values_;
}

Vec6& BeamStressResultants::values() {
    return values_;
}

} // namespace fem
