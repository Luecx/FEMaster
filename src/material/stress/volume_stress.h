#pragma once

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/coordinate_system.h"

namespace fem {

struct VolumeStress {
    enum class Component : int {
        XX = 0,
        YY = 1,
        ZZ = 2,
        YZ = 3,
        XZ = 4,
        XY = 5
    };

    VolumeStress() = default;
    explicit VolumeStress(const Vec6& voigt);
    explicit VolumeStress(const Mat3& tensor);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();
    [[nodiscard]] Mat3        tensor() const;

    [[nodiscard]] VolumeStress transformed(const cos::Basis& from_basis,
                                            const cos::Basis& to_basis) const;

    static Mat6 get_transformation_matrix(const cos::Basis& from_basis,
                                          const cos::Basis& to_basis);

protected:
    Vec6 voigt_{Vec6::Zero()};
};

} // namespace fem
