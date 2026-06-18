/**
 * @file strain.h
 * @brief Declares the strain tensor utility in Voigt form.
 *
 * Provides transformation helpers to rotate strain vectors between coordinate
 * systems.
 *
 * @see src/material/strain.cpp
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
 * @struct Strain
 * @brief Represents a symmetric strain tensor in engineering Voigt notation.
 *
 * Stored Voigt ordering:
 *
 *     [e11, e22, e33, gamma23, gamma13, gamma12]^T
 */
struct Strain {
    enum class Component : int {
        XX = 0,
        YY = 1,
        ZZ = 2,
        GammaYZ = 3,
        GammaZX = 4,
        GammaXY = 5
    };

    enum class ComponentDerived : int {
        TrueYZ,
        TrueZX,
        TrueXY
    };

    Strain() = default;
    explicit Strain(const Vec6& voigt);
    explicit Strain(const Mat3& tensor);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;
    Precision  operator[](ComponentDerived component) const;

    void set(ComponentDerived component, Precision value);

    [[nodiscard]] const Vec6& voigt() const;
    [[nodiscard]] Vec6&       voigt();
    [[nodiscard]] Mat3        tensor() const;

    [[nodiscard]] Strain transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;

    /**
     * @brief Computes the 6x6 strain transformation matrix between two bases.
     */
    static Mat6 get_transformation_matrix(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    );

    /**
     * @brief Legacy helper for a basis matrix interpreted as the full transform.
     */
    static Mat6 get_transformation_matrix(const cos::Basis& basis);

protected:
    Vec6 voigt_{Vec6::Zero()};
};

struct LinearizedStrain : Strain {
    LinearizedStrain() = default;
    explicit LinearizedStrain(const Vec6& voigt);
    explicit LinearizedStrain(const Mat3& tensor);

    [[nodiscard]] LinearizedStrain transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;
};

struct GreenLagrangeStrain : Strain {
    GreenLagrangeStrain() = default;
    explicit GreenLagrangeStrain(const Vec6& voigt);
    explicit GreenLagrangeStrain(const Mat3& tensor);

    [[nodiscard]] GreenLagrangeStrain transformed(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;

    static GreenLagrangeStrain from_deformation_gradient(const Mat3& deformation_gradient);
};

using Strains = std::vector<Strain>; ///< Collection alias for strain vectors.
} // namespace fem
