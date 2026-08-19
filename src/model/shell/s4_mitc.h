/**
 * @file s4_mitc.h
 * @brief Declares the four-node MITC4 shell element.
 *
 * The element uses bilinear quadrilateral geometry together with the classical
 * MITC4 edge-midpoint interpolation of the two transverse-shear components.
 */

#pragma once

#include "shell_simple.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

struct MITC4
    : DefaultShellElement<4, Surface4,
                          math::quadrature::Domain::DOMAIN_ISO_QUAD,
                          math::quadrature::Order::ORDER_CUBIC>
{
    using Base = DefaultShellElement<4, Surface4,
                                     math::quadrature::Domain::DOMAIN_ISO_QUAD,
                                     math::quadrature::Order::ORDER_CUBIC>;

    using LocalCoords     = typename Base::LocalCoords;
    using ShapeFunction   = typename Base::ShapeFunction;
    using ShapeDerivative = typename Base::ShapeDerivative;
    using Jacobian        = typename Base::Jacobian;

    MITC4(ID id, std::array<ID, 4> nodes)
        : Base(id, nodes) {}

    // Preserve the MITC4 dynamic type while cloning only the persistent shell
    // topology. All compiled assembly state is attached afterwards.
    ElementPtr copy() const override { return std::make_shared<MITC4>(elem_id, node_ids); }

    std::string type_name() const override { return "MITC4"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
        );
    }

    // MITC4: shear strains are tied at edge midpoints and interpolated into
    // the current natural position (r,s).
    StaticMatrix<2, 12>
    strain_disp_shear_at(Precision r, Precision s, const LocalCoords& xy) override
    {
        // Build one shear-strain row at a tying point. Shell transverse shear
        // follows generalized order gamma_xz, gamma_yz.
        auto Bs_at = [&](Precision rr, Precision ss, bool xz) -> StaticMatrix<1, 12> {
            ShapeDerivative dH_rs = this->shape_derivative(rr, ss);
            Jacobian        J     = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy));
            Mat2            invJ  = J.inverse();

            // Transform natural derivatives into the local element xy plane.
            auto dH_xy = (dH_rs * invJ).transpose();
            auto H     = this->shape_function(rr, ss);

            StaticMatrix<1, 12> row;
            row.setZero();

            for (int i = 0; i < 4; ++i) {
                const Index c_w  = 3 * i + 0;
                const Index c_rx = 3 * i + 1;
                const Index c_ry = 3 * i + 2;

                if (xz) {
                    // gamma_xz = dw/dx + theta_y
                    row(0, c_w ) += dH_xy(0, i);
                    row(0, c_ry) += H(i);
                } else {
                    // gamma_yz = dw/dy - theta_x
                    row(0, c_w ) += dH_xy(1, i);
                    row(0, c_rx) += -H(i);
                }
            }
            return row;
        };

        // Edge-midpoint interpolation weights in natural coordinates.
        const Precision a_m = Precision(0.5) * (Precision(1.0) - s);
        const Precision a_p = Precision(0.5) * (Precision(1.0) + s);
        const Precision b_m = Precision(0.5) * (Precision(1.0) - r);
        const Precision b_p = Precision(0.5) * (Precision(1.0) + r);

        // Interpolate the tied shear strains to the requested point.
        StaticMatrix<1, 12> row_xz =
            a_m * Bs_at(Precision(0.0), Precision(-1.0), true) +
            a_p * Bs_at(Precision(0.0), Precision(+1.0), true);

        StaticMatrix<1, 12> row_yz =
            b_m * Bs_at(Precision(-1.0), Precision(0.0), false) +
            b_p * Bs_at(Precision(+1.0), Precision(0.0), false);

        StaticMatrix<2, 12> Bs;
        Bs.setZero();
        Bs.row(0) = row_xz;
        Bs.row(1) = row_yz;
        return Bs;
    }
};

} // namespace fem::model
