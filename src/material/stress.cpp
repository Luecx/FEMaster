/******************************************************************************
 * @file stress.cpp
 * @brief Implements stress rotation utilities in Voigt notation.
 *
 * @see src/material/stress.h
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "stress.h"

namespace fem {

/******************************************************************************
 * @copydoc Stress::transform
 ******************************************************************************/
Stress Stress::transform(const cos::Basis& basis) const {
    return Stress{Stress::get_transformation_matrix(basis) * (*this)};
}

/******************************************************************************
 * @copydoc Stress::get_transformation_matrix
 ******************************************************************************/
Mat6 Stress::get_transformation_matrix(const cos::Basis& basis) {
    Precision R11 = basis(0, 0), R12 = basis(0, 1), R13 = basis(0, 2);
    Precision R21 = basis(1, 0), R22 = basis(1, 1), R23 = basis(1, 2);
    Precision R31 = basis(2, 0), R32 = basis(2, 1), R33 = basis(2, 2);

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

} // namespace fem
