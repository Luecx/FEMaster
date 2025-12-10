#ifndef S4_MITC_H
#define S4_MITC_H

#include "shell_simple.h"

#include "../geometry/surface/surface4.h"
namespace fem::model {
struct MITC4 : DefaultShellElement<4, Surface4,
                      quadrature::Domain::DOMAIN_ISO_QUAD,
                      quadrature::Order::ORDER_CUBIC> // dein 2x2 ist damit ok
{
    using Base = DefaultShellElement<4, Surface4,
                                     quadrature::Domain::DOMAIN_ISO_QUAD,
                                     quadrature::Order::ORDER_CUBIC>;

    MITC4(ID id, std::array<ID,4> nodes) : Base(id, nodes) {}

    std::string type_name() const override { return "MITC4"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID,4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID,4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]});
    }

    // Nur diese Methode ist MITC-spezifisch:
    StaticMatrix<2, 12>
    strain_disp_shear_at(Precision r, Precision s, const LocalCoords& xy) override
    {
        // Tying-Points:
        struct RS { Precision r, s; };
        const RS tp_xz[2] = { {0.0, -1.0}, {0.0, +1.0} };
        const RS tp_yz[2] = { {-1.0, 0.0}, {+1.0, 0.0} };

        auto Bs_at = [&](Precision rr, Precision ss, bool xz)->StaticMatrix<1,12> {
            ShapeDerivative dH_rs = this->shape_derivative(rr, ss);
            Jacobian        J    = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy));
            Mat2 invJ = J.inverse();
            auto dH_xy = (dH_rs * invJ).transpose();  // (2×N)
            auto H     = this->shape_function(rr, ss);

            StaticMatrix<1,12> row; row.setZero();
            for (int i = 0; i < 4; ++i) {
                const Index c_w  = 3*i + 0;
                const Index c_rx = 3*i + 1;
                const Index c_ry = 3*i + 2;
                if (xz) {
                    row(0, c_w ) += dH_xy(0, i); // dw/dx
                    row(0, c_rx) += H(i);        // + rx
                } else {
                    row(0, c_w ) += dH_xy(1, i); // dw/dy
                    row(0, c_ry) += H(i);        // + ry
                }
            }
            return row;
        };

        // Gewichte von Tying-Points -> (r,s)
        const Precision a_m = 0.5 * (1.0 - s); // s=-1
        const Precision a_p = 0.5 * (1.0 + s); // s=+1
        const Precision b_m = 0.5 * (1.0 - r); // r=-1
        const Precision b_p = 0.5 * (1.0 + r); // r=+1

        StaticMatrix<1,12> row_xz = a_m * Bs_at(0.0, -1.0, /*xz=*/true)
                                  + a_p * Bs_at(0.0, +1.0, /*xz=*/true);

        StaticMatrix<1,12> row_yz = b_m * Bs_at(-1.0, 0.0, /*xz=*/false)
                                  + b_p * Bs_at(+1.0, 0.0, /*xz=*/false);

        StaticMatrix<2,12> Bs; Bs.setZero();
        Bs.row(0) = row_xz;
        Bs.row(1) = row_yz;
        return Bs;
    }
};
}

#endif //S4_MITC_H
