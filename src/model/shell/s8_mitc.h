#ifndef S8_MITC_H
#define S8_MITC_H

#include "shell_simple.h"
#include "../geometry/surface/surface8.h"

namespace fem::model {

struct S8MITC
    : DefaultShellElement<
          8,
          Surface8,
          quadrature::Domain::DOMAIN_ISO_QUAD,
          quadrature::Order::ORDER_QUINTIC
      > {
    using Base = DefaultShellElement<
        8,
        Surface8,
        quadrature::Domain::DOMAIN_ISO_QUAD,
        quadrature::Order::ORDER_QUINTIC
    >;

    using LocalCoords     = typename Base::LocalCoords;
    using ShapeFunction   = typename Base::ShapeFunction;
    using ShapeDerivative = typename Base::ShapeDerivative;
    using Jacobian        = typename Base::Jacobian;

    using ShearMatrix = StaticMatrix<2, 24>;
    using ShearRow    = StaticMatrix<1, 24>;

    S8MITC(ID p_elem_id, std::array<ID, 8> p_node)
        : Base(p_elem_id, p_node) {}

    std::string type_name() const override {
        return "S8MITC";
    }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface8>(
            surface_id == 1
                ? std::array<ID, 8>{
                      this->nodes()[0],
                      this->nodes()[1],
                      this->nodes()[2],
                      this->nodes()[3],
                      this->nodes()[4],
                      this->nodes()[5],
                      this->nodes()[6],
                      this->nodes()[7]
                  }
                : std::array<ID, 8>{
                      this->nodes()[0],
                      this->nodes()[3],
                      this->nodes()[2],
                      this->nodes()[1],
                      this->nodes()[7],
                      this->nodes()[6],
                      this->nodes()[5],
                      this->nodes()[4]
                  }
        );
    }

    ShearMatrix strain_disp_shear_at(
        Precision          r,
        Precision          s,
        const LocalCoords& xy
    ) override {
        auto lagrange_minus = [](Precision q) {
            return Precision(0.5) * q * (q - Precision(1));
        };

        auto lagrange_zero = [](Precision q) {
            return Precision(1) - q * q;
        };

        auto lagrange_plus = [](Precision q) {
            return Precision(0.5) * q * (q + Precision(1));
        };

        auto physical_shear_matrix = [&](
            Precision rr,
            Precision ss
        ) -> ShearMatrix {
            ShapeDerivative dH_rs = this->shape_derivative(rr, ss);
            Jacobian        J     = this->jacobian(
                dH_rs,
                const_cast<LocalCoords&>(xy)
            );

            const Precision detJ = J.determinant();

            logging::error(
                detJ > Precision(0),
                "S8MITC element ",
                this->elem_id,
                " has a non-positive Jacobian determinant at (",
                rr,
                ", ",
                ss,
                "): ",
                detJ
            );

            const auto dH_xy = (dH_rs * J.inverse()).transpose();
            const auto H     = this->shape_function(rr, ss);

            ShearMatrix B_xy;
            B_xy.setZero();

            for (Index i = 0; i < 8; ++i) {
                const Index c_w  = 3 * i;
                const Index c_rx = c_w + 1;
                const Index c_ry = c_w + 2;

                B_xy(0, c_w ) =  dH_xy(0, i);
                B_xy(0, c_ry) =  H(i);

                B_xy(1, c_w ) =  dH_xy(1, i);
                B_xy(1, c_rx) = -H(i);
            }

            return B_xy;
        };

        auto covariant_shear_row = [&](
            Precision rr,
            Precision ss,
            Index     component
        ) -> ShearRow {
            ShapeDerivative dH_rs = this->shape_derivative(rr, ss);
            Jacobian        J     = this->jacobian(
                dH_rs,
                const_cast<LocalCoords&>(xy)
            );

            const ShearMatrix B_xy = physical_shear_matrix(rr, ss);
            const ShearMatrix B_rs = J.transpose() * B_xy;

            return B_rs.row(component);
        };

        const ShearRow B_rz =
              lagrange_minus(s)
            * covariant_shear_row(
                  Precision(0),
                  Precision(-1),
                  0
              )
            + lagrange_zero(s)
            * covariant_shear_row(
                  Precision(0),
                  Precision(0),
                  0
              )
            + lagrange_plus(s)
            * covariant_shear_row(
                  Precision(0),
                  Precision(+1),
                  0
              );

        const ShearRow B_sz =
              lagrange_minus(r)
            * covariant_shear_row(
                  Precision(-1),
                  Precision(0),
                  1
              )
            + lagrange_zero(r)
            * covariant_shear_row(
                  Precision(0),
                  Precision(0),
                  1
              )
            + lagrange_plus(r)
            * covariant_shear_row(
                  Precision(+1),
                  Precision(0),
                  1
              );

        ShearMatrix B_rs;
        B_rs.row(0) = B_rz;
        B_rs.row(1) = B_sz;

        ShapeDerivative dH_eval = this->shape_derivative(r, s);
        Jacobian        J_eval  = this->jacobian(
            dH_eval,
            const_cast<LocalCoords&>(xy)
        );

        const Precision detJ_eval = J_eval.determinant();

        logging::error(
            detJ_eval > Precision(0),
            "S8MITC element ",
            this->elem_id,
            " has a non-positive Jacobian determinant at (",
            r,
            ", ",
            s,
            "): ",
            detJ_eval
        );

        const ShearMatrix B_xy =
            J_eval.inverse().transpose() * B_rs;

        ShearMatrix B;
        B.row(0) = B_xy.row(1);
        B.row(1) = B_xy.row(0);

        return B;
    }
};

} // namespace fem::model

#endif // S8_MITC_H
