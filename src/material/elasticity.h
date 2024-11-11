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

    StaticMatrix<6,6> transformation(StaticMatrix<3,3> R) {
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
        return T_eps;
    }

    StaticMatrix<6,6> transformation_der(StaticMatrix<3,3> R, StaticMatrix<3,3> dR) {
        Precision R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
        Precision R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
        Precision R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);

        Precision dR11 = dR(0, 0), dR12 = dR(0, 1), dR13 = dR(0, 2);
        Precision dR21 = dR(1, 0), dR22 = dR(1, 1), dR23 = dR(1, 2);
        Precision dR31 = dR(2, 0), dR32 = dR(2, 1), dR33 = dR(2, 2);

        StaticMatrix<6, 6> dT_eps;
        dT_eps << 2 * R11 * dR11, 2 * R21 * dR21, 2 * R31 * dR31, R11 * dR21 + dR11 * R21, R21 * dR31 + dR21 * R31, R31 * dR11 + dR31 * R11,
                  2 * R12 * dR12, 2 * R22 * dR22, 2 * R32 * dR32, R12 * dR22 + dR12 * R22, R22 * dR32 + dR22 * R32, R32 * dR12 + dR32 * R12,
                  2 * R13 * dR13, 2 * R23 * dR23, 2 * R33 * dR33, R13 * dR23 + dR13 * R23, R23 * dR33 + dR23 * R33, R33 * dR13 + dR33 * R13,
                  2 * (R11 * dR12 + dR11 * R12), 2 * (R21 * dR22 + dR21 * R22), 2 * (R31 * dR32 + dR31 * R32), dR11 * R22 + R11 * dR22 + dR12 * R21 + R12 * dR21, dR21 * R32 + R21 * dR32 + dR22 * R31 + R22 * dR31, dR31 * R12 + R31 * dR12 + dR32 * R11 + R32 * dR11,
                  2 * (R12 * dR13 + dR12 * R13), 2 * (R22 * dR23 + dR22 * R23), 2 * (R32 * dR33 + dR32 * R33), dR23 * R12 + R23 * dR12 + dR13 * R22 + R13 * dR22, dR22 * R33 + R22 * dR33 + dR32 * R23 + R32 * dR23, dR32 * R13 + R32 * dR13 + dR33 * R12 + R33 * dR12,
                  2 * (R11 * dR13 + dR11 * R13), 2 * (R23 * dR21 + dR23 * R21), 2 * (R33 * dR31 + dR33 * R31), dR13 * R21 + R13 * dR21 + dR11 * R23 + R11 * dR23, dR23 * R31 + R23 * dR31 + dR21 * R33 + R21 * dR33, dR33 * R11 + R33 * dR11 + dR31 * R13 + R31 * dR13;

        return dT_eps;
    }



    template<size_t D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed(StaticMatrix<D, D> R) {
        // Raise a compile-time error if D = 2, as it's not implemented
        static_assert(D == 3, "Only 3D transformation is implemented.");

        // Get elasticity matrix or tensor for dimension D (assume get<D>() returns it)
        auto elasticity = get<D>();
        StaticMatrix<6, 6> T_eps = transformation(R);
        // Perform the final transformation by multiplying with elasticity matrix
        return  T_eps.transpose() * elasticity * T_eps;
    }

    template<size_t D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed_derivative(StaticMatrix<D, D> R, StaticMatrix<D, D> R_der) {
        // Raise a compile-time error if D = 2, as it's not implemented
        static_assert(D == 3, "Only 3D transformation is implemented.");

        StaticMatrix<6,6> T_eps = transformation(R);
        StaticMatrix<6,6> dT_eps = transformation_der(R, R_der);

        return dT_eps.transpose() * get<D>() * T_eps + T_eps.transpose() * get<D>() * dT_eps;
    }
};

using ElasticityPtr = std::shared_ptr<Elasticity>;
}    // namespace _material
}    // namespace fem
