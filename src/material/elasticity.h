/******************************************************************************
 * @file elasticity.h
 * @brief Declares the base interface for material elasticity models.
 *
 * Elasticity models provide stiffness matrices for various configurations and
 * offer transformation utilities for rotated coordinate frames.
 *
 * @see src/material/elasticity.cpp
 * @see src/material/material.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../core/core.h"

#include <memory>

namespace fem {
namespace material {

/******************************************************************************
 * @struct Elasticity
 * @brief Abstract base class describing material stiffness behaviour.
 ******************************************************************************/
struct Elasticity {
    virtual ~Elasticity() = default;

    /// Returns the in-plane (2D) stiffness matrix.
    virtual StaticMatrix<3, 3> get_2d() = 0;

    /// Returns the full 3D stiffness matrix.
    virtual StaticMatrix<6, 6> get_3d() = 0;

    /// Returns the shear stiffness matrix for shell elements.
    virtual StaticMatrix<2, 2> get_shear(Precision thickness) = 0;

    /// Returns the bending stiffness matrix for shell elements.
    virtual StaticMatrix<3, 3> get_bend(Precision thickness) = 0;

    /// Returns the membrane stiffness matrix (defaults to `get_2d`).
    virtual StaticMatrix<3, 3> get_memb() { return get_2d(); }

    /******************************************************************************
     * @brief Retrieves the stiffness matrix for the specified dimension.
     *
     * @tparam D Spatial dimension (2 or 3).
     * @return StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> Stiffness matrix.
     ******************************************************************************/
    template<Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get();

    /******************************************************************************
     * @brief Builds the 3D transformation matrix for a rotation `R`.
     *
     * @param R Rotation matrix.
     * @return StaticMatrix<6, 6> Transformation matrix.
     ******************************************************************************/
    StaticMatrix<6, 6> transformation(StaticMatrix<3, 3> R);

    /******************************************************************************
     * @brief Computes the derivative of the transformation matrix.
     *
     * @param R Rotation matrix.
     * @param dR Derivative of the rotation matrix.
     * @return StaticMatrix<6, 6> Derivative of the transformation matrix.
     ******************************************************************************/
    StaticMatrix<6, 6> transformation_der(StaticMatrix<3, 3> R, StaticMatrix<3, 3> dR);

    /******************************************************************************
     * @brief Transforms the stiffness matrix using a rotation.
     *
     * @tparam D Spatial dimension (must be 3 for transformation).
     * @param R Rotation matrix.
     * @return StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> Transformed stiffness.
     ******************************************************************************/
    template<Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed(StaticMatrix<D, D> R);

    /******************************************************************************
     * @brief Computes the derivative of the transformed stiffness matrix.
     *
     * @tparam D Spatial dimension (must be 3 for transformation).
     * @param R Rotation matrix.
     * @param R_der Derivative of the rotation matrix.
     * @return StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> Derivative of the transformed stiffness.
     ******************************************************************************/
    template<Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed_derivative(StaticMatrix<D, D> R, StaticMatrix<D, D> R_der);

    /******************************************************************************
     * @brief Casts the elasticity instance to a concrete type.
     *
     * @tparam T Target type deriving from `Elasticity`.
     * @return T* Pointer to the instance or `nullptr` if the cast fails.
     ******************************************************************************/
    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
};

using ElasticityPtr = std::shared_ptr<Elasticity>;

} // namespace material
} // namespace fem
