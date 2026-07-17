#include "volume_strain_linearized.h"

namespace fem {

VolumeStrainLinearized::VolumeStrainLinearized(const Vec6& voigt)
    : VolumeStrain(voigt) {}

VolumeStrainLinearized::VolumeStrainLinearized(const Mat3& tensor)
    : VolumeStrain(tensor) {}

VolumeStrainLinearized VolumeStrainLinearized::transformed(const cos::Basis& from_basis,
                                                           const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrainLinearized(transformed);
}

} // namespace fem
