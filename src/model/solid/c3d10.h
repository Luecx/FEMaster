/**
 * @file c3d10.h
 * @brief Declares the ten-node quadratic tetrahedral solid element.
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
 * @brief Ten-node quadratic tetrahedral continuum element.
 *
 * Quadratic tetrahedral interpolation uses four vertex and six edge nodes.
 * Volume and material quadrature are selected separately.
 * The element stores global node identifiers in interpolation order and
 * inherits geometry transformations, constitutive state handling, thermal
 * operators and mechanical evaluation from SolidElement. Concrete methods
 * define shape functions, quadrature and topology-specific recovery. Copies
 * retain connectivity; Model::compile() establishes model and section bindings.
 */
struct C3D10 : public SolidElement<10> {
    // Construction and concrete element identity
    C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds);

    // Recreate only the persistent element topology. Instance-specific ids,
    // section assignment and runtime state are established by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D10>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D10"; }

    // Natural-coordinate interpolation. Node ordering matches connectivity
    // and the columns of the mechanical strain-displacement operator.
    StaticMatrix<10, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<10, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<10, 3> node_coords_local() override;

    // Boundary faces reuse the global connectivity for surface loads and constraints
    SurfacePtr surface(ID surface_id) override;

    // Volume quadrature; a separate stiffness rule defines material-point storage
    const math::quadrature::Quadrature& integration_scheme() const override;
    const math::quadrature::Quadrature& integration_scheme_stiffness() const override;

protected:
    // Constant natural-space recovery operator from constitutive points to nodes
    const RowMatrix& extrapolation_matrix() override {
        static const RowMatrix matrix = math::extrapolate(
            this->stress_strain_ip_rst(), this->node_coords_local(),
            {math::ExtrapolationBasis::F1,
             math::ExtrapolationBasis::FR,
             math::ExtrapolationBasis::FS,
             math::ExtrapolationBasis::FT});
        return matrix;
    }
};
} }
