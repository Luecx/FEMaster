//
// Created by f_eggers on 24.01.2025.
//

#ifndef STRAIN_H
#define STRAIN_H

#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"

namespace fem {

struct Strain : Vec6 {

    [[nodiscard]] Strain transform(const fem::cos::Basis& R) const {
        // Assume this is in Voigt notation
        return Strain{Strain::get_transformation_matrix(R) * (*this)};
    }

    static Mat6 get_transformation_matrix(const fem::cos::Basis& R) {
        // Extract rotation matrix components
        Precision R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
        Precision R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
        Precision R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);

        // Transformation matrix for strain in Voigt notation
        StaticMatrix<6, 6> T_eps;
        T_eps <<
            R11 * R11, R12 * R12, R13 * R13, R12 * R13, R11 * R13, R11 * R12,
            R21 * R21, R22 * R22, R23 * R23, R22 * R23, R21 * R23, R21 * R22,
            R31 * R31, R32 * R32, R33 * R33, R32 * R33, R31 * R33, R31 * R32,
            2 * R21 * R31, 2 * R22 * R32, 2 * R23 * R33, R22 * R33 + R23 * R32, R21 * R33 + R23 * R31, R21 * R32 + R22 * R31,
            2 * R11 * R31, 2 * R12 * R32, 2 * R13 * R33, R12 * R33 + R13 * R32, R11 * R33 + R13 * R31, R11 * R32 + R12 * R31,
            2 * R11 * R21, 2 * R12 * R22, 2 * R13 * R23, R12 * R23 + R13 * R22, R11 * R23 + R13 * R21, R11 * R22 + R12 * R21;

        return T_eps;
    }

};


using Strains = std::vector<Strain>;

}

#endif //STRAIN_H
