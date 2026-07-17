#include "volume_strain_green_lagrange.h"

namespace fem {

VolumeStrainGreenLagrange::VolumeStrainGreenLagrange(const Vec6& voigt)
    : VolumeStrain(voigt) {}

VolumeStrainGreenLagrange::VolumeStrainGreenLagrange(const Mat3& tensor)
    : VolumeStrain(tensor) {}

VolumeStrainGreenLagrange VolumeStrainGreenLagrange::transformed(const cos::Basis& from_basis,
                                                                 const cos::Basis& to_basis) const {
    const Vec6 transformed = get_transformation_matrix(from_basis, to_basis) * voigt_;
    return VolumeStrainGreenLagrange(transformed);
}

VolumeStrainGreenLagrange VolumeStrainGreenLagrange::from_deformation_gradient(
    const Mat3& deformation_gradient
) {
    const Mat3 strain_tensor = Precision(0.5)
        * (deformation_gradient.transpose() * deformation_gradient - Mat3::Identity());
    return VolumeStrainGreenLagrange(strain_tensor);
}

} // namespace fem
