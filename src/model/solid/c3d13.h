/**
 * @file c3d13.h
 * @brief Declares the thirteen-node quadratic pyramid solid element.
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D13 : public SolidElement<13>{
    C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids);

    // Recreate only the persistent element topology. Dense assembly ids and
    // runtime bindings are assigned after cloning by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D13>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D13"; }

    const math::quadrature::Quadrature& integration_scheme() const override;

    SurfacePtr surface(ID surface_id) override;

    StaticMatrix<13, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> node_coords_local() override;

protected:
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
