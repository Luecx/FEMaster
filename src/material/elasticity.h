#pragma once

#include "../core/core.h"
#include <memory>

namespace fem {
namespace material {
struct Elasticity {
    public:
    virtual StaticMatrix<3, 3> get_2d() = 0;
    virtual StaticMatrix<6, 6> get_3d() = 0;

    virtual ~Elasticity() = default;

    template<size_t D>
    StaticMatrix<D == 3 ? 6:3,D == 3 ? 6:3> get() {
        if constexpr (D == 2){
            return get_2d();
        }
        if constexpr (D == 3){
            return get_3d();
        }
        throw std::out_of_range("indexing _material elasticity outside of valid dimensions which is 2 or 3");
    }

    template<size_t D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed(StaticMatrix<D, D> R) {
        // Raise a compile-time error if D = 2, as it's not implemented
        static_assert(D == 3, "Only 3D transformation is implemented.");

        // Get elasticity matrix or tensor for dimension D (assume get<D>() returns it)
        auto               elasticity = get<D>();

        Precision R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
        Precision R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
        Precision R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);

        StaticMatrix<6, 6> T_eps;
        T_eps << R11*R11, R21*R21, R31*R31, R11*R21, R21*R31, R31*R11,
                 R12*R12, R22*R22, R32*R32, R12*R22, R22*R32, R32*R12,
                 R13*R13, R23*R23, R33*R33, R13*R23, R23*R33, R33*R13,
                 2*R11*R12, 2*R21*R22, 2*R31*R32, R11*R22 + R12*R21, R21*R32 + R31*R22, R31*R12 + R32*R11,
                 2*R12*R13, 2*R22*R23, 2*R33*R32, R23*R12 + R13*R22, R22*R33 + R32*R23, R32*R13 + R12*R33,
                 2*R11*R13, 2*R23*R21, 2*R33*R31, R13*R21 + R11*R23, R23*R31 + R21*R33, R33*R11 + R31*R13;

        // Perform the final transformation by multiplying with elasticity matrix
        return  T_eps.transpose() * elasticity * T_eps;
    }
};

using ElasticityPtr = std::shared_ptr<Elasticity>;
}    // namespace _material
}    // namespace fem
