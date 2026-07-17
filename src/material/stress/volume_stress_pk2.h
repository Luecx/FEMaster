#pragma once

#include "volume_stress.h"
#include "volume_stress_cauchy.h"

namespace fem {

struct VolumeStressPK2 : VolumeStress {
    VolumeStressPK2() = default;
    explicit VolumeStressPK2(const Vec6& voigt);
    explicit VolumeStressPK2(const Mat3& tensor);

    [[nodiscard]] VolumeStressPK2 transformed(const cos::Basis& from_basis,
                                              const cos::Basis& to_basis) const;

    [[nodiscard]] VolumeStressCauchy to_cauchy(const Mat3& deformation_gradient) const;
};

} // namespace fem
