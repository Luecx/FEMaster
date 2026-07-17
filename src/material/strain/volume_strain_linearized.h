#pragma once

#include "volume_strain.h"

namespace fem {

struct VolumeStrainLinearized : VolumeStrain {
    VolumeStrainLinearized() = default;
    explicit VolumeStrainLinearized(const Vec6& voigt);
    explicit VolumeStrainLinearized(const Mat3& tensor);

    [[nodiscard]] VolumeStrainLinearized transformed(const cos::Basis& from_basis,
                                                      const cos::Basis& to_basis) const;
};

} // namespace fem
