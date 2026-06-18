/**
* @file stress.h
 * @brief Declares the stress tensor utility in Voigt form.
 *
 * Provides transformation helpers to rotate stress vectors between coordinate
 * systems.
 *
 * @see src/material/stress.cpp
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"

#include <vector>

namespace fem {

/**
 * @struct Stress
 * @brief Represents a symmetric stress tensor in solver Voigt notation.
 *
 * Stored Voigt ordering:
 *
 *     [s11, s22, s33, s23, s13, s12]^T
 */
struct Stress {
    enum class Component : int {
        XX = 0,
        YY = 1,
        ZZ = 2,
        YZ = 3,
        ZX = 4,
        XY = 5
    };

    Stress() = default;
    explicit Stress(const Vec6& voigt);
    explicit Stress(const Mat3& tensor);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();
    [[nodiscard]] Mat3        tensor() const;

    [[nodiscard]] Stress transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;

    /**
     * @brief Computes the 6x6 stress transformation matrix between two bases.
     */
    static Mat6 get_transformation_matrix(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    );

protected:
    Vec6 voigt_{Vec6::Zero()};
};

struct CauchyStress : Stress {
    CauchyStress() = default;
    explicit CauchyStress(const Vec6& voigt);
    explicit CauchyStress(const Mat3& tensor);

    [[nodiscard]] CauchyStress transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;
};

struct PK2Stress : Stress {
    PK2Stress() = default;
    explicit PK2Stress(const Vec6& voigt);
    explicit PK2Stress(const Mat3& tensor);

    [[nodiscard]] PK2Stress transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;

    [[nodiscard]] CauchyStress to_cauchy(const Mat3& deformation_gradient) const;
};

using Stresses = std::vector<Stress>; ///< Collection alias for stress vectors.

} // namespace fem
