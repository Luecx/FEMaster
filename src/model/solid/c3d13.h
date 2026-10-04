/**
 * @file c3d13.h
 * @brief Declares the thirteen-node quadratic pyramid solid element.
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
 * @brief Thirteen-node quadratic pyramid continuum element.
 *
 * Rational interpolation maps the natural pyramid, including its collapsed
 * apex. Boundary surface extraction is currently unavailable.
 * The element stores global node identifiers in interpolation order and
 * inherits geometry transformations, constitutive state handling, thermal
 * operators and mechanical evaluation from SolidElement. Concrete methods
 * define shape functions, quadrature and topology-specific recovery. Copies
 * retain connectivity; Model::compile() establishes model and section bindings.
 */
struct C3D13 : public SolidElement<13> {
    // Construction and concrete element identity
    C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids);

    // Recreate only the persistent element topology. Dense assembly ids and
    // runtime bindings are assigned after cloning by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D13>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D13"; }

    // Volume quadrature; a separate stiffness rule defines material-point storage
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Boundary faces reuse the global connectivity for surface loads and constraints
    SurfacePtr surface(ID surface_id) override;

    // Natural-coordinate interpolation. Node ordering matches connectivity
    // and the columns of the mechanical strain-displacement operator.
    StaticMatrix<13, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> node_coords_local() override;

protected:
    // Constant natural-space recovery operator from constitutive points to nodes
    const RowMatrix& extrapolation_matrix() override {
        static const RowMatrix matrix = math::extrapolate(
            this->stress_strain_ip_rst(), this->node_coords_local(),
            {math::ExtrapolationBasis::F1,
             math::ExtrapolationBasis::FR,
             math::ExtrapolationBasis::FS,
             math::ExtrapolationBasis::FT,
             math::ExtrapolationBasis::FRR,
             math::ExtrapolationBasis::FSS,
             math::ExtrapolationBasis::FTT,
             math::ExtrapolationBasis::FRS,
             math::ExtrapolationBasis::FRT,
             math::ExtrapolationBasis::FST});
        return matrix;
    }
};
} }
