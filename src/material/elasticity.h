/******************************************************************************
* @file elasticity.h
* @brief Defines the Elasticity base class for material elasticity properties.
*
* The Elasticity class is an abstract base class providing the interface for
* calculating stiffness matrices for materials under different conditions
* (e.g., 2D, 3D, shear, and bending). It includes transformation utilities
* for rotating stiffness matrices in 3D.
*
* @see elasticity.cpp
* @note This class uses template methods for dimension-specific operations.
* @author Finn Eggers
* @date 20.12.2024
******************************************************************************/

#pragma once  // Ensures this file is included only once during compilation

#include "../core/core.h"        // Core FEM types
#include <memory>               // For std::shared_ptr

namespace fem::material {

/******************************************************************************
* @class Elasticity
* @brief Abstract base class for material elasticity properties.
*
* Provides an interface for calculating material stiffness matrices in 2D, 3D,
* shear, and bending contexts. Includes utilities for transforming elasticity
* matrices in 3D using rotation matrices.
******************************************************************************/
struct Elasticity {
public:
    virtual ~Elasticity() = default;

    /******************************************************************************
    * @brief Computes the 2D stiffness matrix.
    * @return A 3x3 static matrix representing the 2D stiffness.
    ******************************************************************************/
    virtual StaticMatrix<3, 3> get_2d() = 0;

    /******************************************************************************
    * @brief Computes the 3D stiffness matrix.
    * @return A 6x6 static matrix representing the 3D stiffness.
    ******************************************************************************/
    virtual StaticMatrix<6, 6> get_3d() = 0;

    /******************************************************************************
    * @brief Computes the shear stiffness matrix for shell elements.
    * @param t Thickness of the shell element.
    * @return A 2x2 static matrix representing shear stiffness.
    ******************************************************************************/
    virtual StaticMatrix<2, 2> get_shear(Precision t) = 0;

    /******************************************************************************
    * @brief Computes the bending stiffness matrix for shell elements.
    * @param t Thickness of the shell element.
    * @return A 3x3 static matrix representing bending stiffness.
    ******************************************************************************/
    virtual StaticMatrix<3, 3> get_bend(Precision t) = 0;

    /******************************************************************************
    * @brief Computes the membrane stiffness matrix for shell elements.
    * Defaults to the 2D stiffness matrix.
    * @return A 3x3 static matrix representing membrane stiffness.
    ******************************************************************************/
    virtual StaticMatrix<3, 3> get_memb() { return get_2d(); }

    /******************************************************************************
    * @brief Retrieves the stiffness matrix for a given dimension.
    * @tparam D The spatial dimension (2 or 3).
    * @return The stiffness matrix for the given dimension.
    ******************************************************************************/
    template <Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get();

    /******************************************************************************
    * @brief Transforms the stiffness matrix using a rotation matrix in 3D.
    * @param R The 3x3 rotation matrix.
    * @return The transformed 6x6 stiffness matrix.
    ******************************************************************************/
    StaticMatrix<6, 6> transformation(StaticMatrix<3, 3> R);

    /******************************************************************************
    * @brief Computes the derivative of the transformation matrix with respect to rotation.
    * @param R The 3x3 rotation matrix.
    * @param dR The derivative of the 3x3 rotation matrix.
    * @return The derivative of the 6x6 transformation matrix.
    ******************************************************************************/
    StaticMatrix<6, 6> transformation_der(StaticMatrix<3, 3> R, StaticMatrix<3, 3> dR);

    /******************************************************************************
    * @brief Computes the transformed stiffness matrix for a given rotation in 3D.
    * @tparam D The spatial dimension (must be 3).
    * @param R The rotation matrix.
    * @return The transformed stiffness matrix.
    ******************************************************************************/
    template <Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed(StaticMatrix<D, D> R);

    /******************************************************************************
    * @brief Computes the derivative of the transformed stiffness matrix.
    * @tparam D The spatial dimension (must be 3).
    * @param R The rotation matrix.
    * @param R_der The derivative of the rotation matrix.
    * @return The derivative of the transformed stiffness matrix.
    ******************************************************************************/
    template <Dim D>
    StaticMatrix<(D == 3 ? 6 : 3), (D == 3 ? 6 : 3)> get_transformed_derivative(StaticMatrix<D, D> R, StaticMatrix<D, D> R_der);

    /******************************************************************************
    * @brief Dynamically casts this instance to a specific type.
    * @tparam T The target type.
    * @return A pointer to the instance cast to the target type.
    ******************************************************************************/
    template <typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
};

using ElasticityPtr = std::shared_ptr<Elasticity>;

}    // namespace fem::material
