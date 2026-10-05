/**
 * @file c3d6.h
 * @brief Declares the six-node linear wedge solid element.
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
 * @brief Six-node linear wedge continuum element.
 *
 * The natural triangular prism combines linear triangular interpolation with linear interpolation
 * through its thickness.
 * The element stores global node identifiers in interpolation order and
 * inherits geometry transformations, constitutive state handling, thermal
 * operators and mechanical evaluation from SolidElement. Concrete methods
 * define shape functions, quadrature and topology-specific recovery. Copies
 * retain connectivity; Model::compile() establishes model and section bindings.
 */
struct C3D6 : public SolidElement<6> {
    // Construction and concrete element identity
    C3D6(ID p_elem_id, const std::array<ID, 6>& p_node_ids);

    // Recreate only the persistent element topology. Instance-specific ids and
    // runtime bindings are assigned by Model::compile() afterwards.
    ElementPtr copy() const override { return std::make_shared<C3D6>(elem_id, node_ids); }

    // Volume quadrature; a separate stiffness rule defines material-point storage
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Boundary faces reuse the global connectivity for surface loads and constraints
    SurfacePtr surface(ID surface_id) override;

    // Natural-coordinate interpolation. Node ordering matches connectivity
    // and the columns of the mechanical strain-displacement operator.
    StaticMatrix<6, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> node_coords_local() override;

    std::string type_name() const override { return "C3D6"; }

protected:
    // Constant natural-space recovery operator from constitutive points to nodes
    const RowMatrix& extrapolation_matrix() override;
};
} }
