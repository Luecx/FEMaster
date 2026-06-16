#pragma once

#include "shell.h"

#include "../geometry/surface/surface4.h"

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>

namespace fem::model {

struct MITC4NL : ShellElement<4> {
    using NodeCoords      = StaticMatrix<4, 3>;
    using LocalCoords     = StaticMatrix<4, 2>;
    using ShapeFunction   = StaticVector<4>;
    using ShapeDerivative = StaticMatrix<4, 2>;
    using Matrix24        = StaticMatrix<24, 24>;
    using Vector24        = StaticVector<24>;
    using Vector8         = StaticVector<8>;
    using Matrix6x24      = StaticMatrix<6, 24>;
    using Matrix8x24      = StaticMatrix<8, 24>;
    using Matrix3x24      = StaticMatrix<3, 24>;
    using Matrix2x24      = StaticMatrix<2, 24>;

    Surface4               geometry;
    quadrature::Quadrature integration_scheme_;

    MITC4NL(ID elem_id, std::array<ID, 4> node_ids)
        : ShellElement<4>(elem_id, node_ids)
        , geometry(node_ids)
        , integration_scheme_(quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_CUBIC) {}

    std::string type_name() const override { return "MITC4NL"; }

    SurfacePtr surface(ID surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]});
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    static ShapeFunction shape_function(Precision r, Precision s) {
        ShapeFunction N;
        N << Precision(0.25) * (Precision(1) - r) * (Precision(1) - s),
             Precision(0.25) * (Precision(1) + r) * (Precision(1) - s),
             Precision(0.25) * (Precision(1) + r) * (Precision(1) + s),
             Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
        return N;
    }

    static ShapeDerivative shape_derivative(Precision r, Precision s) {
        ShapeDerivative dN;
        dN << -Precision(0.25) * (Precision(1) - s), -Precision(0.25) * (Precision(1) - r),
               Precision(0.25) * (Precision(1) - s), -Precision(0.25) * (Precision(1) + r),
               Precision(0.25) * (Precision(1) + s),  Precision(0.25) * (Precision(1) + r),
              -Precision(0.25) * (Precision(1) + s),  Precision(0.25) * (Precision(1) - r);
        return dN;
    }

    static StaticMatrix<4, 3> nodal_rst() {
        StaticMatrix<4, 3> rst;
        rst << -Precision(1), -Precision(1), Precision(0),
                Precision(1), -Precision(1), Precision(0),
                Precision(1),  Precision(1), Precision(0),
               -Precision(1),  Precision(1), Precision(0);
        return rst;
    }

    static Vec3 row_vec(const NodeCoords& coords, Index row) {
        return coords.row(row).transpose();
    }

    Mat3 element_basis(const NodeCoords& coords) {
        constexpr Precision tol = Precision(1e-12);

        const Vec3 x1 = row_vec(coords, 0);
        const Vec3 x2 = row_vec(coords, 1);
        const Vec3 x3 = row_vec(coords, 2);
        const Vec3 x4 = row_vec(coords, 3);

        Vec3 e1 = x2 - x1;
        logging::error(e1.norm() > tol, "MITC4NL: degenerate edge 1-2 in element ", this->elem_id);
        e1.normalize();

        Vec3 e3 = (x3 - x1).cross(x4 - x2);
        if (e3.norm() <= tol) {
            e3 = (x2 - x1).cross(x4 - x1);
        }
        logging::error(e3.norm() > tol, "MITC4NL: degenerate shell normal in element ", this->elem_id);
        e3.normalize();

        Vec3 e2 = e3.cross(e1);
        logging::error(e2.norm() > tol, "MITC4NL: invalid local basis in element ", this->elem_id);
        e2.normalize();

        Mat3 R;
        R.col(0) = e1;
        R.col(1) = e2;
        R.col(2) = e3;
        return R;
    }

    LocalCoords local_coords(const NodeCoords& coords, const Mat3& basis) {
        LocalCoords xy;
        const Vec3  x0 = row_vec(coords, 0);
        for (Index i = 0; i < 4; ++i) {
            const Vec3 local = basis.transpose() * (row_vec(coords, i) - x0);
            xy(i, 0) = local(0);
            xy(i, 1) = local(1);
        }
        return xy;
    }

    static Mat2 jacobian(const ShapeDerivative& dN, const LocalCoords& xy) {
        Mat2 J;
        J(0, 0) = dN.col(0).dot(xy.col(0));
        J(1, 0) = dN.col(0).dot(xy.col(1));
        J(0, 1) = dN.col(1).dot(xy.col(0));
        J(1, 1) = dN.col(1).dot(xy.col(1));
        return J;
    }

    static StaticMatrix<2, 4> shape_derivative_xy(Precision r, Precision s, const LocalCoords& xy) {
        const ShapeDerivative dN = shape_derivative(r, s);
        const Mat2            J  = jacobian(dN, xy);
        return (dN * J.inverse()).transpose();
    }

    static Precision area_weight(Precision r, Precision s, const LocalCoords& xy) {
        return std::abs(jacobian(shape_derivative(r, s), xy).determinant());
    }

    Vec3 global_position(Precision r, Precision s) {
        const ShapeFunction N      = shape_function(r, s);
        const NodeCoords    coords = this->node_coords_global();
        Vec3                x      = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * row_vec(coords, i);
        }
        return x;
    }

    Matrix24 global_to_local_transform(const Mat3& basis) {
        Matrix24 T = Matrix24::Zero();
        for (Index i = 0; i < 4; ++i) {
            T.template block<3, 3>(6 * i,     6 * i    ) = basis.transpose();
            T.template block<3, 3>(6 * i + 3, 6 * i + 3) = basis.transpose();
        }
        return T;
    }

    Vector24 local_dofs(const Field& displacement, const Mat3& basis) {
        Vector24 u = Vector24::Zero();
        for (Index i = 0; i < 4; ++i) {
            const Vec6 row = displacement.row_vec6(static_cast<Index>(this->node_ids[i]));
            u.template segment<3>(6 * i)     = basis.transpose() * row.template head<3>();
            u.template segment<3>(6 * i + 3) = basis.transpose() * row.template tail<3>();
        }
        return u;
    }

    NodeCoords reference_coords_from_current(const Field& displacement) {
        NodeCoords coords = this->node_coords_global();
        for (Index i = 0; i < 4; ++i) {
            const Vec3 u = displacement.row_vec6(static_cast<Index>(this->node_ids[i])).template head<3>();
            coords.row(i) -= u.transpose();
        }
        return coords;
    }

    static Matrix3x24 membrane_B(const StaticMatrix<2, 4>& dNxy) {
        Matrix3x24 B = Matrix3x24::Zero();
        for (Index i = 0; i < 4; ++i) {
            B(0, 6 * i    ) = dNxy(0, i);
            B(1, 6 * i + 1) = dNxy(1, i);
            B(2, 6 * i    ) = dNxy(1, i);
            B(2, 6 * i + 1) = dNxy(0, i);
        }
        return B;
    }

    static Matrix3x24 bending_B(const StaticMatrix<2, 4>& dNxy) {
        Matrix3x24 B = Matrix3x24::Zero();
        for (Index i = 0; i < 4; ++i) {
            B(0, 6 * i + 4) = -dNxy(0, i);
            B(1, 6 * i + 3) =  dNxy(1, i);
            B(2, 6 * i + 3) =  dNxy(0, i);
            B(2, 6 * i + 4) = -dNxy(1, i);
        }
        return B;
    }

    static StaticMatrix<1, 24> shear_row_at(Precision r, Precision s, const LocalCoords& xy, bool xz) {
        const ShapeFunction      N    = shape_function(r, s);
        const StaticMatrix<2, 4> dNxy = shape_derivative_xy(r, s, xy);
        StaticMatrix<1, 24>     row  = StaticMatrix<1, 24>::Zero();

        for (Index i = 0; i < 4; ++i) {
            if (xz) {
                row(0, 6 * i + 2) += dNxy(0, i);
                row(0, 6 * i + 4) += N(i);
            } else {
                row(0, 6 * i + 2) += dNxy(1, i);
                row(0, 6 * i + 3) -= N(i);
            }
        }
        return row;
    }

    static Matrix2x24 shear_B_mitc(Precision r, Precision s, const LocalCoords& xy) {
        const Precision wm = Precision(0.5) * (Precision(1) - s);
        const Precision wp = Precision(0.5) * (Precision(1) + s);
        const Precision vm = Precision(0.5) * (Precision(1) - r);
        const Precision vp = Precision(0.5) * (Precision(1) + r);

        Matrix2x24 B = Matrix2x24::Zero();
        B.row(0) = vm * shear_row_at(-Precision(1), Precision(0), xy, false)
                 + vp * shear_row_at( Precision(1), Precision(0), xy, false);
        B.row(1) = wm * shear_row_at(Precision(0), -Precision(1), xy, true)
                 + wp * shear_row_at(Precision(0),  Precision(1), xy, true);
        return B;
    }

    static Vec3 rotated_director(const Vec3& rotation) {
        const Precision angle = rotation.norm();
        const Vec3      e3    = Vec3::UnitZ();
        if (angle <= Precision(1e-12)) {
            return (e3 + rotation.cross(e3)).normalized();
        }
        return (Eigen::AngleAxis<Precision>(angle, rotation / angle) * e3).normalized();
    }

    Vec3 director_at(Precision r, Precision s, const Vector24& u) {
        const ShapeFunction N = shape_function(r, s);
        Vec3                d = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            d += N(i) * rotated_director(u.template segment<3>(6 * i + 3));
        }
        return d.normalized();
    }

    Vec3 director_derivative(Precision r, Precision s, const LocalCoords& xy, const Vector24& u, Index dir) {
        const StaticMatrix<2, 4> dNxy = shape_derivative_xy(r, s, xy);
        Vec3                    dd   = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            dd += dNxy(dir, i) * rotated_director(u.template segment<3>(6 * i + 3));
        }
        return dd;
    }

    Vec3 displacement_derivative(Precision r, Precision s, const LocalCoords& xy, const Vector24& u, Index dir) {
        const StaticMatrix<2, 4> dNxy = shape_derivative_xy(r, s, xy);
        Vec3                    du   = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            du += dNxy(dir, i) * u.template segment<3>(6 * i);
        }
        return du;
    }

    Precision shear_component_nonlinear(Precision r, Precision s, const LocalCoords& xy, const Vector24& u, bool xz) {
        const Vec3 du = displacement_derivative(r, s, xy, u, xz ? 0 : 1);
        const Vec3 g  = (xz ? Vec3::UnitX() : Vec3::UnitY()) + du;
        const Vec3 d  = director_at(r, s, u);
        return g.dot(d);
    }

    Vec2 shear_strain_mitc_nonlinear(Precision r, Precision s, const LocalCoords& xy, const Vector24& u) {
        const Precision wm = Precision(0.5) * (Precision(1) - s);
        const Precision wp = Precision(0.5) * (Precision(1) + s);
        const Precision vm = Precision(0.5) * (Precision(1) - r);
        const Precision vp = Precision(0.5) * (Precision(1) + r);

        Vec2 gamma;
        gamma(0) = vm * shear_component_nonlinear(-Precision(1), Precision(0), xy, u, false)
                 + vp * shear_component_nonlinear( Precision(1), Precision(0), xy, u, false);
        gamma(1) = wm * shear_component_nonlinear(Precision(0), -Precision(1), xy, u, true)
                 + wp * shear_component_nonlinear(Precision(0),  Precision(1), xy, u, true);
        return gamma;
    }

    static Vec2 shear_yz_xz_to_xz_yz(const Vec2& shear) {
        Vec2 reordered;
        reordered << shear(1), shear(0);
        return reordered;
    }

    static Mat2 shear_matrix_yz_xz_to_xz_yz(const Mat2& shear) {
        Mat2 reordered;
        reordered << shear(1, 1), shear(1, 0),
                     shear(0, 1), shear(0, 0);
        return reordered;
    }

    Vector8 generalized_strain8(Precision r, Precision s, const LocalCoords& xy, const Vector24& u, bool nonlinear) {
        const StaticMatrix<2, 4> dNxy = shape_derivative_xy(r, s, xy);

        Vec3 eps          = membrane_B(dNxy) * u;
        Vec3 kappa        = bending_B(dNxy)  * u;
        Vec2 shear_yz_xz  = shear_B_mitc(r, s, xy) * u;

        if (nonlinear) {
            const Vec3 ux = displacement_derivative(r, s, xy, u, 0);
            const Vec3 uy = displacement_derivative(r, s, xy, u, 1);
            const Vec3 gx = Vec3::UnitX() + ux;
            const Vec3 gy = Vec3::UnitY() + uy;

            eps(0) = Precision(0.5) * (gx.dot(gx) - Precision(1));
            eps(1) = Precision(0.5) * (gy.dot(gy) - Precision(1));
            eps(2) = gx.dot(gy);

            const Vec3 dx = director_derivative(r, s, xy, u, 0);
            const Vec3 dy = director_derivative(r, s, xy, u, 1);

            kappa(0) = -dx.dot(gx);
            kappa(1) = -dy.dot(gy);
            kappa(2) = -(dx.dot(gy) + dy.dot(gx));

            shear_yz_xz = shear_strain_mitc_nonlinear(r, s, xy, u);
        }

        Vector8 strain;
        strain.template segment<3>(0) = eps;
        strain.template segment<3>(3) = kappa;
        strain.template segment<2>(6) = shear_yz_xz_to_xz_yz(shear_yz_xz);
        return strain;
    }

    Vec6 generalized_strain(Precision r, Precision s, const LocalCoords& xy, const Vector24& u, bool nonlinear) {
        const Vector8 strain = generalized_strain8(r, s, xy, u, nonlinear);
        Vec6          result = strain.template segment<6>(0);
        return result;
    }

    Precision topo_scale() {
        Precision scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
        return scale;
    }

    StaticMatrix<6, 6> abd_matrix() {
        return topo_scale() * this->get_section()->get_abd();
    }

    Mat2 shear_matrix() {
        return topo_scale() * this->get_section()->get_shear();
    }

    Mat2 shear_matrix_state_order() {
        return shear_matrix_yz_xz_to_xz_yz(shear_matrix());
    }

    Vector8 generalized_resultant8(const Vector8& strain) {
        const Vec6 nm_strain = strain.template segment<6>(0);
        const Vec2 q_strain  = strain.template segment<2>(6);

        Vector8 result;
        result.template segment<6>(0) = abd_matrix() * nm_strain;
        result.template segment<2>(6) = shear_matrix_state_order() * q_strain;
        return result;
    }

    Vector8 generalized_resultant8_from_state(const Field& ip_stress, Index row) {
        logging::error(ip_stress.components >= 8,
                       "MITC4NL: nonlinear stress state requires at least 8 components");

        Vector8 result;
        for (Index i = 0; i < 8; ++i) {
            result(i) = ip_stress(row, i);
        }
        return result;
    }

    Matrix8x24 strain_displacement_tangent(Precision r,
                                           Precision s,
                                           const LocalCoords& xy,
                                           const Vector24& u,
                                           bool nonlinear) {
        if (!nonlinear) {
            const auto      dNxy = shape_derivative_xy(r, s, xy);
            const Matrix2x24 Bs  = shear_B_mitc(r, s, xy);
            Matrix8x24      B    = Matrix8x24::Zero();
            B.template block<3, 24>(0, 0) = membrane_B(dNxy);
            B.template block<3, 24>(3, 0) = bending_B(dNxy);
            B.row(6) = Bs.row(1);
            B.row(7) = Bs.row(0);
            return B;
        }

        Matrix8x24 B = Matrix8x24::Zero();
        for (Index j = 0; j < 24; ++j) {
            const Precision h = Precision(1e-7) * std::max(Precision(1), std::abs(u(j)));
            Vector24        up = u;
            Vector24        um = u;
            up(j) += h;
            um(j) -= h;

            const Vector8 ep = generalized_strain8(r, s, xy, up, true);
            const Vector8 em = generalized_strain8(r, s, xy, um, true);
            B.col(j) = (ep - em) / (Precision(2) * h);
        }
        return B;
    }

    Vector24 internal_force_local_from_state(const LocalCoords& xy,
                                             const Vector24& u,
                                             const Field& ip_stress,
                                             int ip_offset,
                                             bool nonlinear) {
        Vector24 fint = Vector24::Zero();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto      p      = integration_scheme_.get_point(ip);
            const Index     row    = static_cast<Index>(ip_offset) + ip;
            const Precision dA     = p.w * area_weight(p.r, p.s, xy);
            const Matrix8x24 B     = strain_displacement_tangent(p.r, p.s, xy, u, nonlinear);
            const Vector8    state = generalized_resultant8_from_state(ip_stress, row);
            fint += B.transpose() * state * dA;
        }
        return fint;
    }

    Vector24 internal_force_local(const LocalCoords& xy, const Vector24& u, bool nonlinear) {
        Vector24 fint = Vector24::Zero();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto      p      = integration_scheme_.get_point(ip);
            const Precision dA     = p.w * area_weight(p.r, p.s, xy);
            const Vector8    state = generalized_resultant8(generalized_strain8(p.r, p.s, xy, u, nonlinear));
            const Matrix8x24 B     = strain_displacement_tangent(p.r, p.s, xy, u, nonlinear);
            fint += B.transpose() * state * dA;
        }
        return fint;
    }

    Matrix24 consistent_tangent_local(const LocalCoords& xy, const Vector24& u) {
        Matrix24 K = Matrix24::Zero();
        for (Index j = 0; j < 24; ++j) {
            const Precision h = Precision(1e-6) * std::max(Precision(1), std::abs(u(j)));
            Vector24        up = u;
            Vector24        um = u;
            up(j) += h;
            um(j) -= h;

            const Vector24 fp = internal_force_local(xy, up, true);
            const Vector24 fm = internal_force_local(xy, um, true);
            K.col(j) = (fp - fm) / (Precision(2) * h);
        }

        const Precision drill = std::max(K.diagonal().cwiseAbs().maxCoeff(), Precision(1)) * Precision(1e-9);
        for (Index i = 0; i < 4; ++i) {
            K(6 * i + 5, 6 * i + 5) += drill;
        }
        return K;
    }

    Matrix24 stiffness_local(const LocalCoords& xy) {
        Matrix24 K     = Matrix24::Zero();
        const auto ABD = abd_matrix();
        const Mat2 Ks  = shear_matrix();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto      p    = integration_scheme_.get_point(ip);
            const Precision dA   = p.w * area_weight(p.r, p.s, xy);
            const auto      dNxy = shape_derivative_xy(p.r, p.s, xy);

            Matrix6x24 Bg = Matrix6x24::Zero();
            Bg.template block<3, 24>(0, 0) = membrane_B(dNxy);
            Bg.template block<3, 24>(3, 0) = bending_B(dNxy);

            const Matrix2x24 Bs = shear_B_mitc(p.r, p.s, xy);
            K += Bg.transpose() * ABD * Bg * dA;
            K += Bs.transpose() * Ks  * Bs * dA;
        }

        const Precision drill = std::max(K.diagonal().cwiseAbs().maxCoeff(), Precision(1)) * Precision(1e-9);
        for (Index i = 0; i < 4; ++i) {
            K(6 * i + 5, 6 * i + 5) += drill;
        }
        return Precision(0.5) * (K + K.transpose());
    }

    Matrix24 mass_local(const LocalCoords& xy) {
        Matrix24 M           = Matrix24::Zero();
        const Precision rho  = this->get_density(false);
        const Precision h    = this->get_section()->thickness_;
        const Precision im   = rho * h * h * h / Precision(12);

        if (rho <= Precision(0)) {
            return M;
        }

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto          p     = integration_scheme_.get_point(ip);
            const ShapeFunction N     = shape_function(p.r, p.s);
            const Precision     scale = p.w * area_weight(p.r, p.s, xy);

            for (Index i = 0; i < 4; ++i) {
                for (Index j = 0; j < 4; ++j) {
                    const Precision mt = rho * h  * scale * N(i) * N(j);
                    const Precision mr = im       * scale * N(i) * N(j);
                    for (Index d = 0; d < 3; ++d) {
                        M(6 * i + d,     6 * j + d    ) += mt;
                        M(6 * i + 3 + d, 6 * j + 3 + d) += mr;
                    }
                }
            }
        }
        return M;
    }

    MapMatrix stiffness(Precision* buffer) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Matrix24    T     = global_to_local_transform(basis);

        MapMatrix mapped(buffer, 24, 24);
        mapped = T.transpose() * stiffness_local(xy) * T;
        return mapped;
    }

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Matrix24    T     = global_to_local_transform(basis);
        Matrix24          Kg    = Matrix24::Zero();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto      p     = integration_scheme_.get_point(ip);
            const Index     row   = static_cast<Index>(ip_start_idx) + ip;
            const Precision dA    = p.w * area_weight(p.r, p.s, xy);
            const auto      dNxy  = shape_derivative_xy(p.r, p.s, xy);
            Mat2            Nmat;
            Nmat << ip_stress(row, 0), ip_stress(row, 2),
                    ip_stress(row, 2), ip_stress(row, 1);

            for (Index a = 0; a < 4; ++a) {
                const Vec2 ga(dNxy(0, a), dNxy(1, a));
                for (Index b = 0; b < 4; ++b) {
                    const Vec2      gb = Vec2(dNxy(0, b), dNxy(1, b));
                    const Precision k  = ga.dot(Nmat * gb) * dA;
                    for (Index d = 0; d < 3; ++d) {
                        Kg(6 * a + d, 6 * b + d) += k;
                    }
                }
            }
        }

        MapMatrix mapped(buffer, 24, 24);
        mapped = T.transpose() * Kg * T;
        return mapped;
    }

    MapMatrix consistent_tangent(Precision* buffer, const Field& displacement) {
        const NodeCoords coords = reference_coords_from_current(displacement);
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Matrix24    T     = global_to_local_transform(basis);
        const Vector24    u     = local_dofs(displacement, basis);

        MapMatrix mapped(buffer, 24, 24);
        mapped = T.transpose() * consistent_tangent_local(xy, u) * T;
        return mapped;
    }

    MapMatrix mass(Precision* buffer) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Matrix24    T     = global_to_local_transform(basis);

        MapMatrix mapped(buffer, 24, 24);
        mapped = T.transpose() * mass_local(xy) * T;
        return mapped;
    }

    Precision volume() override {
        return integrate_scalar_field(false, [](const Vec3&) { return Precision(1); });
    }

    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Precision   rho   = scale_by_density ? this->get_density(true) : Precision(1);
        const Precision   h     = this->get_section()->thickness_;
        Precision         value = Precision(0);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto          p  = integration_scheme_.get_point(ip);
            const ShapeFunction N  = shape_function(p.r, p.s);
            Vec3                x  = Vec3::Zero();
            for (Index i = 0; i < 4; ++i) {
                x += N(i) * row_vec(coords, i);
            }
            value += field(x) * rho * h * p.w * area_weight(p.r, p.s, xy);
        }
        return value;
    }

    Vec3 integrate_vector_field(bool scale_by_density, const VecField& field) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Precision   rho   = scale_by_density ? this->get_density(true) : Precision(1);
        const Precision   h     = this->get_section()->thickness_;
        Vec3              value = Vec3::Zero();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto          p = integration_scheme_.get_point(ip);
            const ShapeFunction N = shape_function(p.r, p.s);
            Vec3                x = Vec3::Zero();
            for (Index i = 0; i < 4; ++i) {
                x += N(i) * row_vec(coords, i);
            }
            value += field(x) * rho * h * p.w * area_weight(p.r, p.s, xy);
        }
        return value;
    }

    void integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Precision   rho   = scale_by_density ? this->get_density(true) : Precision(1);
        const Precision   h     = this->get_section()->thickness_;

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto          p = integration_scheme_.get_point(ip);
            const ShapeFunction N = shape_function(p.r, p.s);
            Vec3                x = Vec3::Zero();
            for (Index i = 0; i < 4; ++i) {
                x += N(i) * row_vec(coords, i);
            }

            const Vec3 f = field(x) * rho * h * p.w * area_weight(p.r, p.s, xy);
            for (Index i = 0; i < 4; ++i) {
                const Index node = static_cast<Index>(this->node_ids[i]);
                node_loads(node, 0) += N(i) * f(0);
                node_loads(node, 1) += N(i) * f(1);
                node_loads(node, 2) += N(i) * f(2);
            }
        }
    }

    Mat3 integrate_tensor_field(bool scale_by_density, const TenField& field) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Precision   rho   = scale_by_density ? this->get_density(true) : Precision(1);
        const Precision   h     = this->get_section()->thickness_;
        Mat3              value = Mat3::Zero();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto          p = integration_scheme_.get_point(ip);
            const ShapeFunction N = shape_function(p.r, p.s);
            Vec3                x = Vec3::Zero();
            for (Index i = 0; i < 4; ++i) {
                x += N(i) * row_vec(coords, i);
            }
            value += field(x) * rho * h * p.w * area_weight(p.r, p.s, xy);
        }
        return value;
    }

    RowMatrix stress_strain_nodal_rst() override {
        const auto coords = nodal_rst();
        RowMatrix  rst(4, 3);
        for (Index i = 0; i < 4; ++i) {
            rst.row(i) = coords.row(i);
        }
        return rst;
    }

    RowMatrix stress_strain_ip_rst() override {
        RowMatrix rst(integration_scheme_.count(), 3);
        for (Index i = 0; i < integration_scheme_.count(); ++i) {
            rst(i, 0) = integration_scheme_.get_point(i).r;
            rst(i, 1) = integration_scheme_.get_point(i).s;
            rst(i, 2) = Precision(0);
        }
        return rst;
    }

    void stress_strain_at(Precision r,
                          Precision s,
                          Precision zeta,
                          const Field& displacement,
                          bool nonlinear,
                          Vec6& strain_out,
                          Vec6& stress_out) {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Vector24    u     = local_dofs(displacement, basis);
        const Precision   h     = this->get_section()->thickness_;

        const Vector8 qstrain   = generalized_strain8(r, s, xy, u, nonlinear);
        const Vector8 resultant = generalized_resultant8(qstrain);
        const Vec3    eps       = qstrain.template segment<3>(0);
        const Vec3    kappa     = qstrain.template segment<3>(3);
        const Vec2    gamma     = qstrain.template segment<2>(6);
        const Vec3    N         = resultant.template segment<3>(0);
        const Vec3    M         = resultant.template segment<3>(3);
        const Vec2    Q         = resultant.template segment<2>(6);
        const Precision z       = zeta * h / Precision(2);

        const Vec3 plane_strain = eps + z * kappa;
        const Vec3 plane_stress = N / h + z * (Precision(12) / (h * h * h)) * M;
        const Vec2 shear_stress = Q / h;

        strain_out.setZero();
        stress_out.setZero();
        strain_out(0) = plane_strain(0);
        strain_out(1) = plane_strain(1);
        strain_out(3) = gamma(1);
        strain_out(4) = gamma(0);
        strain_out(5) = plane_strain(2);

        stress_out(0) = plane_stress(0);
        stress_out(1) = plane_stress(1);
        stress_out(3) = shear_stress(1);
        stress_out(4) = shear_stress(0);
        stress_out(5) = plane_stress(2);
    }

    void compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        logging::error(strain != nullptr || stress != nullptr,
                       "MITC4NL: compute_stress_strain requires at least one output field");
        logging::error(rst.cols() >= 3,
                       "MITC4NL: stress/strain evaluation coordinates require at least 3 columns");

        for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
            Vec6 strain_value;
            Vec6 stress_value;
            stress_strain_at(rst(i, 0), rst(i, 1), rst(i, 2), displacement, use_green_lagrange_nl, strain_value, stress_value);

            const Index row = static_cast<Index>(offset) + i;
            if (strain) {
                for (Index c = 0; c < strain->components; ++c) {
                    (*strain)(row, c) = c < 6 ? strain_value(c) : Precision(0);
                }
            }
            if (stress) {
                for (Index c = 0; c < stress->components; ++c) {
                    (*stress)(row, c) = c < 6 ? stress_value(c) : Precision(0);
                }
            }
        }
    }

    void compute_stress_state(Field& stress_state,
                              const Field& displacement,
                              int offset,
                              bool use_green_lagrange_nl) override {
        logging::error(stress_state.components >= 8,
                       "MITC4NL: stress state requires at least 8 components");

        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Vector24    u     = local_dofs(displacement, basis);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto    p      = integration_scheme_.get_point(ip);
            const Vector8 result = generalized_resultant8(generalized_strain8(p.r, p.s, xy, u, use_green_lagrange_nl));
            const Index   row    = static_cast<Index>(offset) + ip;

            for (Index c = 0; c < stress_state.components; ++c) {
                stress_state(row, c) = c < 8 ? result(c) : Precision(0);
            }
        }
    }

    void compute_internal_force_nonlinear(Field& node_forces,
                                          const Field& displacement,
                                          const Field& ip_stress,
                                          int ip_offset) override {
        const NodeCoords coords = reference_coords_from_current(displacement);
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Matrix24    T     = global_to_local_transform(basis);
        const Vector24    u     = local_dofs(displacement, basis);
        const Vector24    fint  = internal_force_local_from_state(xy, u, ip_stress, ip_offset, true);

        const Vector24 fglob = T.transpose() * fint;
        for (Index i = 0; i < 4; ++i) {
            const Index node = static_cast<Index>(this->node_ids[i]);
            for (Index d = 0; d < 6; ++d) {
                node_forces(node, d) += fglob(6 * i + d);
            }
        }
    }

    bool compute_shell_section_forces(Field& section_forces,
                                      Field& contribution_count,
                                      const Field& displacement) override {
        const NodeCoords coords = this->node_coords_global();
        const Mat3       basis  = element_basis(coords);
        const LocalCoords xy    = local_coords(coords, basis);
        const Vector24    u     = local_dofs(displacement, basis);
        const auto        rst   = nodal_rst();

        for (Index i = 0; i < 4; ++i) {
            const Precision r   = rst(i, 0);
            const Precision s   = rst(i, 1);
            const Vector8   gen = generalized_strain8(r, s, xy, u, true);
            const Vector8   res = generalized_resultant8(gen);
            const Index     node = static_cast<Index>(this->node_ids[i]);

            section_forces(node, 0) += res(0);
            section_forces(node, 1) += res(1);
            section_forces(node, 2) += res(2);
            section_forces(node, 3) += res(3);
            section_forces(node, 4) += res(4);
            section_forces(node, 5) += res(5);
            section_forces(node, 6) += res(6);
            section_forces(node, 7) += res(7);
            contribution_count(node, 0) += Precision(1);
        }
        return true;
    }
};

} // namespace fem::model
