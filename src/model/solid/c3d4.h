/**
 * @file C3D4.h
 * @brief C3D4.h defines a 4-node tetrahedral element (C3D4) for solid
 * mechanics within a finite element model (FEM). This element computes
 * shape functions, their derivatives, and uses a suitable quadrature scheme
 * for integration.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 27.08.2024
 *
 *
 * The topology declares natural interpolation, quadrature and recovery data.
 * SolidElement owns the common kinematics, constitutive state interface,
 * mechanical evaluation and thermal operators.
 *
 * @see SolidElement
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

/**
 * @brief Four-node linear tetrahedral continuum element.
 *
 * The natural tetrahedron uses constant shape gradients and one material integration point.
 * The element stores global node identifiers in interpolation order and
 * inherits geometry transformations, constitutive state handling, thermal
 * operators and mechanical evaluation from SolidElement. Concrete methods
 * define shape functions, quadrature and topology-specific recovery. Copies
 * retain connectivity; Model::compile() establishes model and section bindings.
 */
struct C3D4 : public SolidElement<4> {
    // Construction and concrete element identity
    C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds);

    // Recreate only the persistent topology. Model::compile() rewires the
    // returned element to the target Instance and initializes runtime state.
    ElementPtr copy() const override { return std::make_shared<C3D4>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D4"; }

    // Natural-coordinate interpolation. Node ordering matches connectivity
    // and the columns of the mechanical strain-displacement operator.
    StaticMatrix<4, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<4, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<4, 3> node_coords_local() override;

    // Volume quadrature; a separate stiffness rule defines material-point storage
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Boundary faces reuse the global connectivity for surface loads and constraints
    SurfacePtr surface(ID surface_id) override;

protected:
    // Constant natural-space recovery operator from constitutive points to nodes
    const RowMatrix& extrapolation_matrix() override {
        static const RowMatrix matrix = math::extrapolate(
            this->stress_strain_ip_rst(), this->node_coords_local(),
            {math::ExtrapolationBasis::F1});
        return matrix;
    }
};
} } // namespace fem::model
