#pragma once

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/coordinate_system.h"

namespace fem {

struct VolumeStrain {
    enum class Component : int {
        XX      = 0,
        YY      = 1,
        ZZ      = 2,
        GammaYZ = 3,
        GammaXZ = 4,
        GammaXY = 5
    };

    enum class TensorComponent : int {
        YZ,
        XZ,
        XY
    };

    VolumeStrain() = default;
    explicit VolumeStrain(const Vec6& voigt);
    explicit VolumeStrain(const Mat3& tensor);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;
    Precision  operator[](TensorComponent component) const;

    void set(TensorComponent component, Precision value);

    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();
    [[nodiscard]] Mat3        tensor() const;

    [[nodiscard]] VolumeStrain transformed(const cos::Basis& from_basis,
                                            const cos::Basis& to_basis) const;

    static Mat6 get_transformation_matrix(const cos::Basis& from_basis,
                                          const cos::Basis& to_basis);

protected:
    Vec6 voigt_{Vec6::Zero()};
};

} // namespace fem
