#pragma once

#include "volume_strain.h"

namespace fem {

struct VolumeStrainGreenLagrange : VolumeStrain {
    VolumeStrainGreenLagrange() = default;
    explicit VolumeStrainGreenLagrange(const Vec6& voigt);
    explicit VolumeStrainGreenLagrange(const Mat3& tensor);

    [[nodiscard]] VolumeStrainGreenLagrange transformed(const cos::Basis& from_basis,
                                                        const cos::Basis& to_basis) const;

    static VolumeStrainGreenLagrange from_deformation_gradient(const Mat3& deformation_gradient);
};

} // namespace fem
