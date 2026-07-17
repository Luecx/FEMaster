#include "volume_stress_cauchy.h"

namespace fem {

VolumeStressCauchy::VolumeStressCauchy(const Vec6& voigt)
    : VolumeStress(voigt) {}

VolumeStressCauchy::VolumeStressCauchy(const Mat3& tensor)
    : VolumeStress(tensor) {}

VolumeStressCauchy VolumeStressCauchy::transformed(const cos::Basis& from_basis,
                                                   const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStressCauchy(transformed);
}

} // namespace fem
