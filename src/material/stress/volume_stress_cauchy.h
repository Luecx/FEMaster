#pragma once

#include "volume_stress.h"

namespace fem {

struct VolumeStressCauchy : VolumeStress {
    VolumeStressCauchy() = default;
    explicit VolumeStressCauchy(const Vec6& voigt);
    explicit VolumeStressCauchy(const Mat3& tensor);

    [[nodiscard]] VolumeStressCauchy transformed(const cos::Basis& from_basis,
                                                 const cos::Basis& to_basis) const;
};

} // namespace fem
