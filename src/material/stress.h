//
// Created by f_eggers on 24.01.2025.
//

#ifndef STRESS_H
#define STRESS_H

#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"

namespace fem {

struct Stress : Vec6 {

    [[nodiscard]] Stress transform(const fem::cos::Basis& R) const {
        // assume this is in voigt notation
        return Stress{Stress::get_transformation_matrix(R) * (*this)};
    }

    static Mat6 get_transformation_matrix(const fem::cos::Basis& R) {
        // Extract rotation matrix components
        Precision R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
        Precision R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
        Precision R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);

        // Transformation matrix for stress in Voigt notation
        StaticMatrix<6, 6> T_sigma;
        T_sigma <<
            R11 * R11, R12 * R12, R13 * R13, 2 * R12 * R13, 2 * R11 * R13, 2 * R11 * R12,
            R21 * R21, R22 * R22, R23 * R23, 2 * R22 * R23, 2 * R21 * R23, 2 * R21 * R22,
            R31 * R31, R32 * R32, R33 * R33, 2 * R32 * R33, 2 * R31 * R33, 2 * R31 * R32,
            R21 * R31, R22 * R32, R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
            R11 * R31, R12 * R32, R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
            R11 * R21, R12 * R22, R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;

        return T_sigma;
    }

};


}

#endif //STRESS_H
