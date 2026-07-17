#pragma once

#include "../../core/types_eig.h"

namespace fem {

struct BeamStressResultants {
    BeamStressResultants() = default;
    explicit BeamStressResultants(const Vec6& values);

    [[nodiscard]] const Vec6& values() const;
    [[nodiscard]] Vec6&       values();

private:
    Vec6 values_{Vec6::Zero()};
};

} // namespace fem
