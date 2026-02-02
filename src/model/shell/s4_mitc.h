#ifndef S4_MITC_H
#define S4_MITC_H

#include "shell_simple.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

struct MITC4
    : DefaultShellElement<4, Surface4,
                          quadrature::Domain::DOMAIN_ISO_QUAD,
                          quadrature::Order::ORDER_CUBIC>
{
    using Base = DefaultShellElement<4, Surface4,
                                     quadrature::Domain::DOMAIN_ISO_QUAD,
                                     quadrature::Order::ORDER_CUBIC>;

    using LocalCoords     = typename Base::LocalCoords;
    using ShapeFunction   = typename Base::ShapeFunction;
    using ShapeDerivative = typename Base::ShapeDerivative;
    using Jacobian        = typename Base::Jacobian;

    MITC4(ID id, std::array<ID, 4> nodes)
        : Base(id, nodes) {}

    std::string type_name() const override { return "MITC4"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
        );
    }

    // MITC4: shear strains are "tied" at edge midpoints and interpolated into (r,s)
    StaticMatrix<2, 12>
    strain_disp_shear_at(Precision r, Precision s, const LocalCoords& xy) override
    {
        // Helper: build 1x12 row for a shear strain component at (rr,ss)
        // In YOUR convention (see DefaultShellElement::strain_disp_shear):
        //   gamma_xz = dw/dx + theta_y   (-> uses ry with + sign)
        //   gamma_yz = dw/dy - theta_x   (-> uses rx with - sign)
        auto Bs_at = [&](Precision rr, Precision ss, bool xz) -> StaticMatrix<1, 12> {
            ShapeDerivative dH_rs = this->shape_derivative(rr, ss);
            Jacobian        J     = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy));
            Mat2            invJ  = J.inverse();

            // dN/dx, dN/dy in local element xy plane
            auto dH_xy = (dH_rs * invJ).transpose();   // (2 x 4): row0=dN/dx, row1=dN/dy
            auto H     = this->shape_function(rr, ss); // (4)

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

        // Tying points (edge midpoints) in (r,s) space:
        //   gamma_xz tied on s = +/-1 at r = 0
        //   gamma_yz tied on r = +/-1 at s = 0
        const Precision a_m = Precision(0.5) * (Precision(1.0) - s); // weight for s=-1
        const Precision a_p = Precision(0.5) * (Precision(1.0) + s); // weight for s=+1

        const Precision b_m = Precision(0.5) * (Precision(1.0) - r); // weight for r=-1
        const Precision b_p = Precision(0.5) * (Precision(1.0) + r); // weight for r=+1

        // Interpolate tied shear strains to (r,s)
        StaticMatrix<1, 12> row_xz =
            a_m * Bs_at(Precision(0.0), Precision(-1.0), /*xz=*/true) +
            a_p * Bs_at(Precision(0.0), Precision(+1.0), /*xz=*/true);

        StaticMatrix<1, 12> row_yz =
            b_m * Bs_at(Precision(-1.0), Precision(0.0), /*xz=*/false) +
            b_p * Bs_at(Precision(+1.0), Precision(0.0), /*xz=*/false);

        StaticMatrix<2, 12> Bs;
        Bs.setZero();
        Bs.row(0) = row_xz;
        Bs.row(1) = row_yz;
        return Bs;
    }
};

} // namespace fem::model

#endif // S4_MITC_H
