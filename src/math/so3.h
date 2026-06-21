/**
 * @file so3.h
 * @brief Utility functions for finite 3D rotations.
 *
 * Provides the exponential-map rotation matrix for a rotation vector and its
 * first and second derivatives with respect to the rotation-vector components.
 */

#pragma once

#include "../core/types_eig.h"

#include <array>

namespace fem::so3 {

Mat3 rotation_matrix(const Vec3& rotation_vector);

void rotation_matrix_first_derivatives(
    const Vec3&          rotation_vector,
    Mat3&                rotation_matrix,
    std::array<Mat3, 3>& first_derivatives
);

void rotation_matrix_second_derivatives(
    const Vec3&                         rotation_vector,
    Mat3&                               rotation_matrix,
    std::array<Mat3, 3>&                first_derivatives,
    std::array<std::array<Mat3, 3>, 3>& second_derivatives
);

} // namespace fem::so3
