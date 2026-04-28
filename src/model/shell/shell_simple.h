//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_SIMPLE_H
#define SHELL_SIMPLE_H

#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"
#include "../geometry/surface/surface.h"
#include "../geometry/surface/surface3.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface6.h"
#include "../geometry/surface/surface8.h"
#include "shell.h"

// this file is partially based on the implementation of
// https://github.com/JWock82/Pynite/blob/main/Archived/MITC4.py
//

namespace fem::model {
template<Index N, typename SFType, quadrature::Domain INT_D, quadrature::Order INT_O>
struct DefaultShellElement : public ShellElement<N> {
    SFType                 geometry;
    quadrature::Quadrature integration_scheme_;

    using Axes            = Mat3;
    using NodeCoords      = StaticMatrix<N, 3>;
    using LocalCoords     = StaticMatrix<N, 2>;
    using Jacobian        = StaticMatrix<2, 2>;
    using ShapeDerivative = StaticMatrix<N, 2>;
    using ShapeFunction   = StaticVector<N>;

    DefaultShellElement(ID p_elem_id, std::array<ID, N> p_node)
        : ShellElement<N>(p_elem_id, p_node)
        , geometry(p_node)
        , integration_scheme_(INT_D, INT_O) {}
    ~DefaultShellElement() override = default;

    // when assuming a planar element, the normal is constant
    // we can also project the shell and use a local x,y system
    // this must not be confused with the local r,s system which is used for integration
    Axes get_xyz_axes() {
        NodeCoords node_coords = this->node_coords_global();

        Vec3       n1          = node_coords.row(0);
        Vec3       n2          = node_coords.row(1);
        Vec3       n3          = node_coords.row(2);

        Vec3       n12         = n2 - n1;
        Vec3       n13         = n3 - n1;

        Vec3       x_axis      = n12 / n12.norm();
        Vec3       z_axis      = n12.cross(n13) / (n12.cross(n13)).norm();
        Vec3       y_axis      = z_axis.cross(x_axis);

        // for each node calculate the local x,y coordinates by projecting onto the axes
        StaticMatrix<3, 3> res;
        res.row(0) = x_axis.transpose();
        res.row(1) = y_axis.transpose();
        res.row(2) = z_axis.transpose();
        return res;
    }

    // Assumption is that the element is planar
    // deviations of that will produce inaccurate results
    // finding the transformation matrix is straight forward
    // the x-axis will be aligned with the direction from node 1 to node 2
    // the z-axis will be perpendicular to the plane of the element
    // the y-axis will be the cross product of the other two
    StaticMatrix<N * 6, N * 6> transformation(Mat3& axes) {
        StaticMatrix<N * 6, N * 6> res = StaticMatrix<N * 6, N * 6>::Zero();
        for (int i = 0; i < 2 * N; ++i) {    // Loop over the 8 nodes
            res.template block<3, 3>(3 * i, 3 * i) = axes;
        }
        return res;
    }

    LocalCoords get_xy_coords(Mat3& axes) {
        NodeCoords node_coords = this->node_coords_global();

        Vec3       x_axis      = axes.row(0).transpose();
        Vec3       y_axis      = axes.row(1).transpose();

        // for each node calculate the local x,y coordinates by projecting onto the axes
        StaticMatrix<N, 2> xy_coords;
        for (Index i = 0; i < N; i++) {
            Vec3 ni         = node_coords.row(i);
            Vec3 n0         = node_coords.row(0);
            Vec3 n          = ni - n0;
            xy_coords(i, 0) = n.dot(x_axis);
            xy_coords(i, 1) = n.dot(y_axis);
        }
        return xy_coords;
    }

    Vec3 global_point(Precision r, Precision s) {
        ShapeFunction H = this->shape_function(r, s);
        NodeCoords node_coords = this->node_coords_global();
        Vec3 point = Vec3::Zero();
        for (Index i = 0; i < N; ++i) {
            point += H(i) * node_coords.row(i).transpose();
        }
        return point;
    }

    Mat3 element_basis_global(Mat3& axes) {
        Mat3 basis;
        basis.col(0) = axes.row(0).transpose();
        basis.col(1) = axes.row(1).transpose();
        basis.col(2) = axes.row(2).transpose();
        return basis;
    }

    Mat3 section_basis_global(Precision r, Precision s, Mat3& axes) {
        auto section = this->get_section();
        if (!section->orientation_) {
            return element_basis_global(axes);
        }

        const Vec3 point_global = global_point(r, s);
        const Vec3 point_local = section->orientation_->to_local(point_global);
        const Mat3 orientation_axes = section->orientation_->get_axes(point_local);

        const Vec3 normal = axes.row(2).transpose().normalized();
        Vec3 n1 = orientation_axes.col(0);
        n1 -= normal * n1.dot(normal);
        if (n1.squaredNorm() < Precision(1e-24)) {
            n1 = orientation_axes.col(1);
            n1 -= normal * n1.dot(normal);
        }
        if (n1.squaredNorm() < Precision(1e-24)) {
            n1 = axes.row(0).transpose();
        }
        n1.normalize();

        Vec3 n2 = normal.cross(n1).normalized();

        Mat3 basis;
        basis.col(0) = n1;
        basis.col(1) = n2;
        basis.col(2) = normal;
        return basis;
    }

    Mat2 section_to_element_2d(Precision r, Precision s, Mat3& axes) {
        const Mat3 section_basis = section_basis_global(r, s, axes);
        const Vec3 ex = axes.row(0).transpose();
        const Vec3 ey = axes.row(1).transpose();

        Mat2 A;
        A(0, 0) = ex.dot(section_basis.col(0));
        A(1, 0) = ey.dot(section_basis.col(0));
        A(0, 1) = ex.dot(section_basis.col(1));
        A(1, 1) = ey.dot(section_basis.col(1));
        return A;
    }

    StaticMatrix<3, 3> plane_strain_transform(const Mat2& section_to_element) {
        const Precision a1 = section_to_element(0, 0);
        const Precision a2 = section_to_element(1, 0);
        const Precision b1 = section_to_element(0, 1);
        const Precision b2 = section_to_element(1, 1);

        StaticMatrix<3, 3> T;
        T << a1 * a1,       a2 * a2,       a1 * a2,
             b1 * b1,       b2 * b2,       b1 * b2,
             2 * a1 * b1,   2 * a2 * b2,   a1 * b2 + a2 * b1;
        return T;
    }

    Mat2 shear_strain_transform(const Mat2& section_to_element) {
        return section_to_element.transpose();
    }

    StaticMatrix<3, 3> transform_plane_stiffness_to_element(const StaticMatrix<3, 3>& section_stiffness,
                                                            const Mat2& section_to_element) {
        const auto T = plane_strain_transform(section_to_element);
        return T.transpose() * section_stiffness * T;
    }

    Mat2 transform_shear_stiffness_to_element(const Mat2& section_stiffness,
                                              const Mat2& section_to_element) {
        const auto T = shear_strain_transform(section_to_element);
        return T.transpose() * section_stiffness * T;
    }

    StaticMatrix<6, 6> transform_abd_stiffness_to_element(const StaticMatrix<6, 6>& section_abd,
                                                          const Mat2& section_to_element) {
        const auto T_plane = plane_strain_transform(section_to_element);
        StaticMatrix<6, 6> T = StaticMatrix<6, 6>::Zero();
        T.template block<3, 3>(0, 0) = T_plane;
        T.template block<3, 3>(3, 3) = T_plane;
        return T.transpose() * section_abd * T;
    }

    ShapeFunction shape_function(Precision r, Precision s) {
        return geometry.shape_function(r, s);
    }

    ShapeDerivative shape_derivative(Precision r, Precision s) {
        return geometry.shape_derivative(r, s);
    }

    // function which relates derivatives in the local r,s to the local x,y system
    // the entries will contain the mapping
    // [dx/dr, dx/ds]
    // [dy/dr, dy/ds]
    // x = N1*x1 + N2*x2 + N3*x3 + N4*x4
    // y = N1*y1 + N2*y2 + N3*y3 + N4*y4
    // hence the derivatives are
    // dx/dr = dN1/dr*x1 + dN2/dr*x2 + dN3/dr*x3 + dN4/dr*x4
    Jacobian jacobian(StaticMatrix<N, 2>& shape_derivatives, StaticMatrix<N, 2>& xy_coords) {
        // computing the jacobian is done by computing the derivatives of the shape functions
        // with the local x,y system
        Jacobian jacobian;
        jacobian(0, 0) = shape_derivatives.col(0).dot(xy_coords.col(0));    // dx/dr
        jacobian(1, 0) = shape_derivatives.col(0).dot(xy_coords.col(1));    // dy/dr

        jacobian(0, 1) = shape_derivatives.col(1).dot(xy_coords.col(0));    // dx/ds
        jacobian(1, 1) = shape_derivatives.col(1).dot(xy_coords.col(1));    // dy/ds
        return jacobian;
    }

    // Integrate a vector field across the shell mid-surface (times thickness)
    void integrate_vec_field(Field& node_loads,
                             bool scale_by_density,
                             const typename ShellElement<N>::VecField& field) override {
        // Local planar frame and projected nodal XY coordinates
        Mat3 axes = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();
        const Precision t = this->get_section()->thickness_;

        Precision rho = 1.0;
        if (scale_by_density) {
            rho = this->get_density(true);
        }

        // For global x at IP, use true global node positions
        auto node_coords_glob = this->node_coords_global();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H   = this->shape_function(r, s); // (N)
            ShapeDerivative dH  = this->shape_derivative(r, s);
            Jacobian        Jrs = this->jacobian(dH, xy_coords); // 2x2
            const Precision detJ = std::abs(Jrs.determinant());

            // Global evaluation point on mid-surface
            Vec3 x_ip = Vec3::Zero();
            for (Index i = 0; i < N; ++i) x_ip += H(i) * node_coords_glob.row(i);

            Vec3 f_ip = field(x_ip) * (rho * t * w * detJ);

            for (Index i = 0; i < N; ++i) {
                const ID n_id = this->node_ids[i];
                const Precision a = H(i);
                node_loads(n_id, 0) += a * f_ip(0);
                node_loads(n_id, 1) += a * f_ip(1);
                node_loads(n_id, 2) += a * f_ip(2);
            }
        }
    }

    // Scalar, vector, and tensor integrations over shell area*thickness
    Precision integrate_scalar_field(bool scale_by_density,
                                     const typename ShellElement<N>::ScalarField& field) override {
        Mat3 axes = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();
        const Precision t = this->get_section()->thickness_;

        Precision rho = 1.0;
        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();
        Precision result = Precision(0);
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H   = this->shape_function(r, s);
            ShapeDerivative dH  = this->shape_derivative(r, s);
            Jacobian        Jrs = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();
            for (Index i = 0; i < N; ++i) x_ip += H(i) * node_coords_glob.row(i);

            result += field(x_ip) * (rho * t * w * detJ);
        }
        return result;
    }

    Vec3 integrate_vector_field(bool scale_by_density,
                                const typename ShellElement<N>::VecField& field) override {
        Mat3 axes = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();
        const Precision t = this->get_section()->thickness_;

        Precision rho = 1.0;
        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();
        Vec3 result = Vec3::Zero();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H   = this->shape_function(r, s);
            ShapeDerivative dH  = this->shape_derivative(r, s);
            Jacobian        Jrs = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();
            for (Index i = 0; i < N; ++i) x_ip += H(i) * node_coords_glob.row(i);

            result += field(x_ip) * (rho * t * w * detJ);
        }
        return result;
    }

    Mat3 integrate_tensor_field(bool scale_by_density,
                                const typename ShellElement<N>::TenField& field) override {
        Mat3 axes = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();
        const Precision t = this->get_section()->thickness_;

        Precision rho = 1.0;
        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();
        Mat3 result = Mat3::Zero();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H   = this->shape_function(r, s);
            ShapeDerivative dH  = this->shape_derivative(r, s);
            Jacobian        Jrs = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();
            for (Index i = 0; i < N; ++i) x_ip += H(i) * node_coords_glob.row(i);

            result += field(x_ip) * (rho * t * w * detJ);
        }
        return result;
    }

    StaticMatrix<3, N * 3> strain_disp_bending(ShapeDerivative& shape_der, Jacobian& jacobian) {
        Mat2                   inv = jacobian.inverse();
        auto                   dH  = (shape_der * inv).transpose();

        StaticMatrix<3, N * 3> res {};

        // dofs are displacement in z, rotation around x, rotation around y
        // its only a function of the rotational dofs (second derivative of displacement)
        for (int i = 0; i < N; i++) {
            res(0, 3 * i)     = 0;
            res(0, 3 * i + 1) = 0;
            res(0, 3 * i + 2) = -dH(0, i);

            res(1, 3 * i)     = 0;
            res(1, 3 * i + 1) = dH(1, i);
            res(1, 3 * i + 2) = 0;

            res(2, 3 * i)     = 0;
            res(2, 3 * i + 1) = dH(0, i);
            res(2, 3 * i + 2) = -dH(1, i);
        }
        // res <<
        //     0,    0,     -dH(0, 0), 0,    0,     -dH(0, 1), 0,    0,     -dH(0, 2), 0,    0,     -dH(0, 3),
        //     0, dH(1, 0),     0,     0, dH(1, 1),     0,     0, dH(1, 2),     0,     0, dH(1, 3),     0    ,
        //     0, dH(0, 0), -dH(1, 0), 0, dH(0, 1), -dH(1, 1), 0, dH(0, 2), -dH(1, 2), 0, dH(0, 3), -dH(1, 3);
        return res;
    }

    virtual StaticMatrix<2, 3 * N>
        strain_disp_shear(ShapeFunction& shape_func, ShapeDerivative& shape_der, Jacobian& jacobian) {
        auto H  = shape_func;
        auto dH = (shape_der * jacobian.inverse()).transpose();

        // StaticMatrix<2, 12> res{};
        // // dofs are displacement in z, rotation around x, rotation around y
        // // its only a function of the rotational dofs (second derivative of displacement)
        // res <<
        //     dH(0,0),  0, H(0), dH(0,1),  0, H(1), dH(0,2),  0, H(2), dH(0,3),  0, H(3),
        //     dH(1,0),-H(0),  0, dH(1,1),-H(1),  0, dH(1,2),-H(2),  0, dH(1,3),-H(3), 0;
        // return res;
        StaticMatrix<2, 3 * N> res {};
        for (int i = 0; i < N; i++) {
            res(0, 3 * i)     = dH(0, i);
            res(0, 3 * i + 1) = 0;
            res(0, 3 * i + 2) = H(i);

            res(1, 3 * i)     = dH(1, i);
            res(1, 3 * i + 1) = -H(i);
            res(1, 3 * i + 2) = 0;
        }
        return res;
    }

    virtual StaticMatrix<2, 3 * N> strain_disp_shear_at(Precision r, Precision s, const LocalCoords& xy_coords) {
        // Default = bisherige (nicht-MITC) Variante
        ShapeFunction   H  = this->shape_function(r, s);
        ShapeDerivative dH = this->shape_derivative(r, s);
        Jacobian        J  = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));
        return this->strain_disp_shear(H, dH, J);
    }

    StaticMatrix<3, 2 * N> strain_disp_membrane(ShapeDerivative& shape_der, Jacobian& jacobian) {
        auto dH = (shape_der * jacobian.inverse()).transpose();
        // StaticMatrix<3, 8> res{};
        // res << dH(0,0), 0      , dH(0,1), 0      , dH(0,2), 0      , dH(0,3), 0      ,
        //            0  , dH(1,0), 0      , dH(1,1), 0      , dH(1,2), 0      , dH(1,3),
        //        dH(1,0), dH(0,0), dH(1,1), dH(0,1), dH(1,2), dH(0,2), dH(1,3), dH(0,3);
        StaticMatrix<3, 2 * N> res {};
        for (int i = 0; i < N; i++) {
            res(0, 2 * i)     = dH(0, i);
            res(0, 2 * i + 1) = 0;

            res(1, 2 * i)     = 0;
            res(1, 2 * i + 1) = dH(1, i);

            res(2, 2 * i)     = dH(1, i);
            res(2, 2 * i + 1) = dH(0, i);
        }
        return res;
    }

    StaticMatrix<6 * N, 6 * N> stiffness_abd_shear(LocalCoords& xy_coords) {
        auto section = this->get_section();
        auto mat_abd_section   = section->get_abd();
        auto mat_shear_section = section->get_shear();
        Mat3 axes = get_xyz_axes();

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
        mat_abd_section *= topo_scale;
        mat_shear_section *= topo_scale;

        std::function<StaticMatrix<5 * N, 5 * N>(Precision, Precision, Precision)> func_abd =
            [this, &mat_abd_section, &xy_coords, &axes](Precision r, Precision s, Precision t) -> StaticMatrix<5 * N, 5 * N> {
            ShapeDerivative shape_der = this->shape_derivative(r, s);
            Jacobian        jac       = this->jacobian(shape_der, xy_coords);
            Precision       jac_det   = jac.determinant();
            auto            Bm        = this->strain_disp_membrane(shape_der, jac);
            auto            Bb        = this->strain_disp_bending(shape_der, jac);
            const Mat2      A         = this->section_to_element_2d(r, s, axes);
            auto            E         = this->transform_abd_stiffness_to_element(mat_abd_section, A);

            StaticMatrix<6, 5 * N> B;
            B.setZero();
            B.template block<3, 2 * N>(0, 0) = Bm;
            B.template block<3, 3 * N>(3, 2 * N) = Bb;
            return B.transpose() * (E * B) * jac_det;
        };

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_shear =
            [this, &mat_shear_section, &xy_coords, &axes](Precision r, Precision s, Precision t) -> StaticMatrix<3 * N, 3 * N> {
            Precision jac_det;                                             // nur fürs Interface hier nicht mehr nötig
            auto      Bs = this->strain_disp_shear_at(r, s, xy_coords);    // <-- NEU
            // wir brauchen dennoch det(J) fürs dA:
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        J    = this->jacobian(dH, xy_coords);
            Precision       detJ = J.determinant();
            const Mat2      A    = this->section_to_element_2d(r, s, axes);
            auto            E    = this->transform_shear_stiffness_to_element(mat_shear_section, A);
            return Bs.transpose() * (E * Bs) * detJ;
        };

        StaticMatrix<5 * N, 5 * N> stiff_abd   = this->integration_scheme().integrate(func_abd);
        StaticMatrix<3 * N, 3 * N> stiff_shear = this->integration_scheme().integrate(func_shear);

        Precision                  k_drill     = 1e32;
        for (int i = 0; i < N; i++) {
            k_drill = std::min(k_drill, stiff_shear(3 * i + 1, 3 * i + 1));
            k_drill = std::min(k_drill, stiff_shear(3 * i + 2, 3 * i + 2));
        }
        k_drill /= 1000;

        StaticMatrix<6 * N, 6 * N> res;
        res.setZero();

        int abd_dof_map[] {
            0,    // x
            1,    // y
            2,    // z
            3,    // rx
            4,    // ry
        };

        for (Index i = 0; i < 5 * N; i++) {
            for (Index j = 0; j < 5 * N; j++) {
                auto i_local_id                 = i < 2 * N ? i / 2 : (i - 2 * N) / 3;
                auto j_local_id                 = j < 2 * N ? j / 2 : (j - 2 * N) / 3;
                auto i_dof                      = i < 2 * N ? i % 2 : 2 + (i - 2 * N) % 3;
                auto j_dof                      = j < 2 * N ? j % 2 : 2 + (j - 2 * N) % 3;

                auto i_glob_index               = i_local_id * 6 + abd_dof_map[i_dof];
                auto j_glob_index               = j_local_id * 6 + abd_dof_map[j_dof];

                res(i_glob_index, j_glob_index) += stiff_abd(i, j);
            }
        }

        int shear_dof_map[] {
            2,    // z
            3,    // rx
            4,    // ry
        };

        for (Index i = 0; i < 3 * N; i++) {
            for (Index j = 0; j < 3 * N; j++) {
                auto i_local_id                 = i / 3;
                auto j_local_id                 = j / 3;
                auto i_dof                      = i % 3;
                auto j_dof                      = j % 3;

                auto i_glob_index               = i_local_id * 6 + shear_dof_map[i_dof];
                auto j_glob_index               = j_local_id * 6 + shear_dof_map[j_dof];

                res(i_glob_index, j_glob_index) += stiff_shear(i, j);
            }
        }

        for (int i = 0; i < N; i++) {
            res(6 * i + 5, 6 * i + 5) = k_drill;
        }

        return res;
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    MapMatrix stiffness(Precision* buffer) override {
        // compute axes and local coordinates
        auto      axes      = get_xyz_axes();
        auto      xy_coords = get_xy_coords(axes);

        MapMatrix mapped {buffer, 6 * N, 6 * N};
        auto      trans = transformation(axes);

        auto stiff = stiffness_abd_shear(xy_coords);
        mapped     = trans.transpose() * stiff * trans;
        return mapped;
    }

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        // 1) lokale Lamina-Achsen & xy-Koordinaten (wie bei Steifigkeit/Masse)
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        // 2) Integrand: ∫ G^T Nmat G detJ dA, nur auf (w)-DOFs im [w,rx,ry]-Block
        Index                                                                      ip_counter = 0;

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_geo =
            [this, &ip_stress, ip_start_idx, &ip_counter, &xy_coords](Precision r,
                                                                      Precision s,
                                                                      Precision /*t*/) -> StaticMatrix<3 * N, 3 * N> {
            // --- Jacobian & dN/dx,y ---
            ShapeDerivative dH_rs = this->shape_derivative(r, s);    // (N×2): dN/dr, dN/ds
            Jacobian        J     = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy_coords));
            Mat2            Jinv  = J.inverse();
            Precision       detJ  = J.determinant();
            auto            dH_xy = (dH_rs * Jinv).transpose();    // (2×N): [dN/dx; dN/dy]

            // --- Membran-Resultierende am IP (aus ip_stress) ---
            // ip_stress: [xx, yy, zz, yz, zx, xy] (wie bei Solids)
            const Index ip_row = static_cast<Index>(ip_start_idx) + ip_counter++;
            const Vec6 v = ip_stress.row_vec6(ip_row);

            const Precision h = this->get_section()->thickness_;
            // Falls ip_stress bereits Resultierende enthält, setze scale=1
            const Precision scale = 1.0;    // oder 1.0, je nachdem was du speicherst
            const Precision nxx   = scale * v(0);
            const Precision nyy   = scale * v(1);
            const Precision nxy   = scale * v(5);

            // Nmat (2x2)
            Eigen::Matrix<Precision, 2, 2> Nmat;
            Nmat << nxx, nxy, nxy, nyy;

            // --- G-Operator (2 × 3N) auf [w,rx,ry]-Blöcke:
            // Zeile 0 mappt ∂w/∂x, Zeile 1 mappt ∂w/∂y; nur die w-Spalten (Index 3*i+0) sind belegt.
            StaticMatrix<2, 3 * N> G;
            G.setZero();
            for (Index i = 0; i < N; ++i) {
                const Index c_w = 3 * i + 0;
                G(0, c_w)       = dH_xy(0, i);    // dw/dx
                G(1, c_w)       = dH_xy(1, i);    // dw/dy
            }

            // --- Lokaler 3N×3N-Block (nur w–w koppelt) ---
            StaticMatrix<3 * N, 3 * N> Kloc = G.transpose() * (Nmat * G) * detJ;
            return Kloc;
        };

        // 3) Integrieren (nur der 3N×3N [w,rx,ry]-Teil)
        StaticMatrix<3 * N, 3 * N> Kg_bend_local = this->integration_scheme().integrate(func_geo);

        // 4) In 6N-Abbildung (z, rx, ry) einsortieren – identisch zu deiner Biege/Schub-Assembly
        StaticMatrix<6 * N, 6 * N> Kg6;
        Kg6.setZero();
        int dof_map[] {2, 3, 4};    // z, rx, ry
        for (Index i = 0; i < 3 * N; ++i) {
            for (Index j = 0; j < 3 * N; ++j) {
                const Index i_node = i / 3, i_ldof = i % 3;
                const Index j_node = j / 3, j_ldof = j % 3;
                const Index I = 6 * i_node + dof_map[i_ldof];
                const Index J = 6 * j_node + dof_map[j_ldof];
                Kg6(I, J)     = Kg_bend_local(i, j);
            }
        }

        // 5) (optional) kein Drill-Anteil in Kg

        // 6) In globale Koordinaten drehen (wie bei 'stiffness')
        MapMatrix mapped {buffer, 6 * N, 6 * N};
        auto      T = transformation(axes);    // 6N×6N
        mapped      = T.transpose() * Kg6 * T;
        return mapped;
    }

    // Hilfsfunktion: ∫_A N N^T det(J) dA  (N×N)
    StaticMatrix<N, N> integrate_NNt(const LocalCoords& xy_coords) {
        std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> fn =
            [this, &xy_coords](Precision r, Precision s, Precision /*t*/) -> StaticMatrix<N, N> {
            ShapeFunction   H    = this->shape_function(r, s);      // (N)
            ShapeDerivative dH   = this->shape_derivative(r, s);    // (N×2), dN/dr, dN/ds
            Jacobian        J    = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));
            Precision       detJ = J.determinant();
            return (H * H.transpose()) * detJ;    // (N×N)
        };
        return this->integration_scheme().integrate(fn);
    }

    MapMatrix mass(Precision* buffer) override {
        // Lokale Achsen & XY-Koordinaten wie bei der Steifigkeit
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        // Material-/Sektionsdaten
        const Precision rho = this->get_density(false);
        if (rho <= Precision(0)) {
            MapMatrix mapped(buffer, 6 * N, 6 * N);
            mapped.setZero();
            return mapped;
        }
        const Precision h   = this->get_section()->thickness_;
        const Precision mt  = rho * h;                     // Translationsmasse pro Fläche
        const Precision mr  = rho * (h * h * h / 12.0);    // Rotations-„Masse“ (Flächenträgheit)

        // Flächenintegral der NN^T-Matrix
        StaticMatrix<N, N> M_NN = integrate_NNt(xy_coords);

        // Lokale (x,y,z, rx,ry,rz) konsistente Massenmatrix (6N × 6N)
        StaticMatrix<6 * N, 6 * N> Mloc;
        Mloc.setZero();

        // 1) Translationsblock: ⊗ I3 auf (ux,uy,uz)
        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mt * M_NN(i, j);
                // (ux,uy,uz) sitzen auf DOFs 0,1,2 im 6er Paket
                for (int k = 0; k < 3; ++k) {
                    const Index I = 6 * i + k;
                    const Index J = 6 * j + k;
                    Mloc(I, J) += mij;
                }
            }
        }

        // 2) Rotationsblock für (rx,ry): ⊗ I2
        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mr * M_NN(i, j);
                // rx -> dof 3, ry -> dof 4
                for (int k = 0; k < 2; ++k) {
                    const Index I = 6 * i + 3 + k;
                    const Index J = 6 * j + 3 + k;
                    Mloc(I, J) += mij;
                }
            }
        }

        // 3) (Optional) sehr kleiner Drill-Massenanteil für rz (dof 5), damit numerisch stabil
        //    Hier nehmen wir z.B. 1e-6 der mittleren Knotentranslationsmasse.
        {
            Precision avg_m_node = 0.0;
            for (Index i = 0; i < N; ++i)
                avg_m_node += mt * M_NN(i, i);
            avg_m_node /= std::max<Precision>(1.0, static_cast<Precision>(N));
            const Precision m_drill = 1e-6 * avg_m_node;
            for (Index i = 0; i < N; ++i) {
                const Index I = 6 * i + 5;    // rz
                Mloc(I, I) += m_drill;
            }
        }

        // In globale Koordinaten drehen (wie bei der Steifigkeit)
        MapMatrix mapped(buffer, 6 * N, 6 * N);
        auto      T = transformation(axes);    // (6N×6N), setzt pro 3er-Block die Achsen
        mapped      = T.transpose() * Mloc * T;
        return mapped;
    }

    Precision volume() override {
        // Dicke aus der Section
        const Precision h = this->get_section()->thickness_;

        // Globale Knotentabelle (POSITION) aus dem Modell
        logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");
        const auto& node_coords_system = *this->_model_data->positions;

        // Flächeninhalt über die Surface-Geometrie (macht 3D-Jacobian×cross-Produkt)
        const Precision A = geometry.area(node_coords_system);

        // Volumen = Dicke * Fläche
        return h * A;
    }

    virtual Stress stress(Field& displacement, Vec3& xyz) {
        Precision r = xyz(0);
        Precision s = xyz(1);
        Precision t = xyz(2);

        Vec6      res;
        res = Vec6::Zero();

        // get the displacement vectors
        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        // get the axes for transformation from global to local
        Mat3 axes = get_xyz_axes();

        for (int i = 0; i < N; i++) {
            ID   node_id             = this->nodes()[i];

            Vec6 displacement_glob   = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz            = displacement_glob.head(3);
            Vec3 disp_rot            = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;   // global → local
            disp_rot = axes * disp_rot;   // global → local

            disp_membrane(2 * i)     = disp_xyz(0);
            disp_membrane(2 * i + 1) = disp_xyz(1);

            disp_bending(3 * i)      = disp_xyz(2);
            disp_bending(3 * i + 1)  = disp_rot(0);
            disp_bending(3 * i + 2)  = disp_rot(1);

            disp_shear(3 * i)        = disp_xyz(2);
            disp_shear(3 * i + 1)    = disp_rot(0);
            disp_shear(3 * i + 2)    = disp_rot(1);
        }

        Precision h = this->get_section()->thickness_;
        auto abd = this->get_section()->get_abd();
        Mat2 mat_shear    = this->get_section()->get_shear();

		// scale material matrices by topo stiffness
        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
		abd          *= topo_scale;
		mat_shear    *= topo_scale;

        // membrane stress
        LocalCoords     xy_coords  = get_xy_coords(axes);
        ShapeFunction   shape_func = this->shape_function(r, s);
        ShapeDerivative shape_der  = this->shape_derivative(r, s);
        Jacobian        jac        = this->jacobian(shape_der, xy_coords);

        // strain displacement matrices
        auto B_membrane = this->strain_disp_membrane(shape_der, jac);
        auto B_bending  = this->strain_disp_bending(shape_der, jac);
        auto B_shear    = this->strain_disp_shear(shape_func, shape_der, jac);

        const Mat2 A = this->section_to_element_2d(r, s, axes);
        const auto T_plane = this->plane_strain_transform(A);
        const auto T_shear = this->shear_strain_transform(A);

        Vec3 eps_membrane = T_plane * (B_membrane * disp_membrane);
        Vec3 kappa        = T_plane * (B_bending  * disp_bending);
        Vec6 generalized_strain;
        generalized_strain << eps_membrane, kappa;
        Vec6 generalized_resultant = abd * generalized_strain;

        const Precision z = t * h / 2;
        const Vec3 membrane_force = generalized_resultant.template segment<3>(0);
        const Vec3 bending_moment = generalized_resultant.template segment<3>(3);
        const Vec3 stress_plane = membrane_force / h + z * (Precision(12) / (h * h * h)) * bending_moment;
        res(0) += stress_plane(0);
        res(1) += stress_plane(1);
        res(5) += stress_plane(2);

        // shear stress
        Vec2 stress_shear = (mat_shear / h) * (T_shear * (B_shear * disp_shear));
        res(3) += stress_shear(0);
        res(4) += stress_shear(1);

        return Stress {res};
    }

    void compute_stress_strain_nodal(Field& displacement,
                                     Field& stress,
                                     Field& strain) override {
        // --- Precompute axes and transformation ---
        Mat3 axes   = get_xyz_axes();
        Mat3 axes_T = axes.transpose(); // local -> global

        // --- Precompute thickness + material matrices ---
        Precision h = this->get_section()->thickness_;

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        auto abd = this->get_section()->get_abd();
        Mat2 mat_shear    = this->get_section()->get_shear();

        abd *= topo_scale;
        mat_shear *= topo_scale;

        // --- Build displacement vectors ONCE (global -> local) ---
        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending.setZero();
        disp_shear.setZero();

        for (int a = 0; a < N; ++a) {
            ID node_id = this->nodes()[a];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            // global -> local only once
            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            disp_membrane(2 * a    ) = disp_xyz(0);
            disp_membrane(2 * a + 1) = disp_xyz(1);

            disp_bending(3 * a    ) = disp_xyz(2);
            disp_bending(3 * a + 1) = disp_rot(0);
            disp_bending(3 * a + 2) = disp_rot(1);

            // current MITC-like shear uses same dofs as bending
            disp_shear(3 * a    ) = disp_xyz(2);
            disp_shear(3 * a + 1) = disp_rot(0);
            disp_shear(3 * a + 2) = disp_rot(1);
        }

        // --- Geometry: local node coordinates & xy coords once ---
        StaticMatrix<N, 2> coords    = this->geometry.node_coords_local();
        LocalCoords        xy_coords = get_xy_coords(axes);

        // --- Loop over nodes and evaluate stress/strain at top/bottom ---
        for (int i = 0; i < N; ++i) {
            ID node_id = this->nodes()[i];

            Precision r = coords(i, 0);
            Precision s = coords(i, 1);

            // shape & jacobian at this (r,s): same for top/bottom
            ShapeFunction   shape_func = this->shape_function(r, s);
            ShapeDerivative shape_der  = this->shape_derivative(r, s);
            Jacobian        jac        = this->jacobian(shape_der, xy_coords);

            auto B_membrane = this->strain_disp_membrane(shape_der, jac);
            auto B_bending  = this->strain_disp_bending(shape_der, jac);
            auto B_shear    = this->strain_disp_shear(shape_func, shape_der, jac);

            const Mat2 A = this->section_to_element_2d(r, s, axes);
            const auto T_plane = this->plane_strain_transform(A);
            const auto T_shear = this->shear_strain_transform(A);

            Vec3 eps_membrane = T_plane * (B_membrane * disp_membrane);
            Vec3 kappa        = T_plane * (B_bending  * disp_bending);
            Vec2 gamma_shear  = T_shear * (B_shear    * disp_shear);

            Vec6 generalized_strain;
            generalized_strain << eps_membrane, kappa;
            Vec6 generalized_resultant = abd * generalized_strain;
            Vec3 membrane_force = generalized_resultant.template segment<3>(0);
            Vec3 bending_moment = generalized_resultant.template segment<3>(3);

            Vec2 stress_shear = (mat_shear / h) * gamma_shear;

            Precision z_top = +0.5 * h;
            Precision z_bot = -0.5 * h;

            Vec6 stress_top_local = Vec6::Zero();
            Vec6 stress_bot_local = Vec6::Zero();
            Vec3 stress_top_plane = membrane_force / h + z_top * (Precision(12) / (h * h * h)) * bending_moment;
            Vec3 stress_bot_plane = membrane_force / h + z_bot * (Precision(12) / (h * h * h)) * bending_moment;

            stress_top_local(0) = stress_top_plane(0);
            stress_top_local(1) = stress_top_plane(1);
            stress_top_local(3) = stress_shear(0);
            stress_top_local(4) = stress_shear(1);
            stress_top_local(5) = stress_top_plane(2);

            stress_bot_local(0) = stress_bot_plane(0);
            stress_bot_local(1) = stress_bot_plane(1);
            stress_bot_local(3) = stress_shear(0);
            stress_bot_local(4) = stress_shear(1);
            stress_bot_local(5) = stress_bot_plane(2);

            // --- Total strain on top and bottom ---
            Vec6 strain_mid_local = Vec6::Zero();
            strain_mid_local(0) = eps_membrane(0);
            strain_mid_local(1) = eps_membrane(1);
            strain_mid_local(3) = gamma_shear(0);
            strain_mid_local(4) = gamma_shear(1);
            strain_mid_local(5) = eps_membrane(2);

            // --- Total strain on top and bottom ---
            Vec6 strain_top_local = strain_mid_local;
            Vec6 strain_bot_local = strain_mid_local;

            strain_top_local(0) += z_top * kappa(0);
            strain_top_local(1) += z_top * kappa(1);
            strain_top_local(5) += z_top * kappa(2);

            strain_bot_local(0) += z_bot * kappa(0);
            strain_bot_local(1) += z_bot * kappa(1);
            strain_bot_local(5) += z_bot * kappa(2);

            // choose side based on stress norm
            Stress stress_top_obj{stress_top_local};
            Stress stress_bot_obj{stress_bot_local};

            bool use_top = stress_top_obj.norm() > stress_bot_obj.norm();

            Vec6 stress_nodal_local = use_top ? stress_top_local : stress_bot_local;
            Vec6 strain_nodal_local = use_top ? strain_top_local : strain_bot_local;

            // stress: transform to global
            Stress stress_nodal_global = Stress{stress_nodal_local}.transform(axes_T);
            Strain strain_nodal_global = Strain{strain_nodal_local}.transform(axes_T);

            const Index node_idx = static_cast<Index>(node_id);
            for (int j = 0; j < 6; ++j) {
                stress(node_idx, j) += stress_nodal_global(j);
                strain(node_idx, j) += strain_nodal_local(j);
            }
        }
    }

    bool supports_shell_stress_surfaces() const override { return true; }

    void compute_shell_stress_surfaces_nodal(Field& displacement,
                                             Field& stress_top,
                                             Field& stress_bot) override {
        Mat3 axes   = get_xyz_axes();
        Mat3 axes_T = axes.transpose();
        Precision h = this->get_section()->thickness_;

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        auto abd = this->get_section()->get_abd();
        Mat2 mat_shear    = this->get_section()->get_shear();

        abd *= topo_scale;
        mat_shear *= topo_scale;

        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending.setZero();
        disp_shear.setZero();

        for (int a = 0; a < N; ++a) {
            ID node_id = this->nodes()[a];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            disp_membrane(2 * a    ) = disp_xyz(0);
            disp_membrane(2 * a + 1) = disp_xyz(1);

            disp_bending(3 * a    ) = disp_xyz(2);
            disp_bending(3 * a + 1) = disp_rot(0);
            disp_bending(3 * a + 2) = disp_rot(1);

            disp_shear(3 * a    ) = disp_xyz(2);
            disp_shear(3 * a + 1) = disp_rot(0);
            disp_shear(3 * a + 2) = disp_rot(1);
        }

        StaticMatrix<N, 2> coords    = this->geometry.node_coords_local();
        LocalCoords        xy_coords = get_xy_coords(axes);

        for (int i = 0; i < N; ++i) {
            ID node_id = this->nodes()[i];

            Precision r = coords(i, 0);
            Precision s = coords(i, 1);

            ShapeFunction   shape_func = this->shape_function(r, s);
            ShapeDerivative shape_der  = this->shape_derivative(r, s);
            Jacobian        jac        = this->jacobian(shape_der, xy_coords);

            auto B_membrane = this->strain_disp_membrane(shape_der, jac);
            auto B_bending  = this->strain_disp_bending(shape_der, jac);
            auto B_shear    = this->strain_disp_shear_at(r, s, xy_coords);

            const Mat2 A = this->section_to_element_2d(r, s, axes);
            const auto T_plane = this->plane_strain_transform(A);
            const auto T_shear = this->shear_strain_transform(A);

            Vec3 eps_membrane = T_plane * (B_membrane * disp_membrane);
            Vec3 kappa        = T_plane * (B_bending  * disp_bending);
            Vec2 gamma_shear  = T_shear * (B_shear    * disp_shear);

            Vec6 generalized_strain;
            generalized_strain << eps_membrane, kappa;
            Vec6 generalized_resultant = abd * generalized_strain;
            Vec3 membrane_force = generalized_resultant.template segment<3>(0);
            Vec3 bending_moment = generalized_resultant.template segment<3>(3);
            Vec2 stress_shear = (mat_shear / h) * gamma_shear;

            Precision z_top = +0.5 * h;
            Precision z_bot = -0.5 * h;

            Vec6 stress_top_local = Vec6::Zero();
            Vec6 stress_bot_local = Vec6::Zero();
            Vec3 stress_top_plane = membrane_force / h + z_top * (Precision(12) / (h * h * h)) * bending_moment;
            Vec3 stress_bot_plane = membrane_force / h + z_bot * (Precision(12) / (h * h * h)) * bending_moment;

            stress_top_local(0) = stress_top_plane(0);
            stress_top_local(1) = stress_top_plane(1);
            stress_top_local(3) = stress_shear(0);
            stress_top_local(4) = stress_shear(1);
            stress_top_local(5) = stress_top_plane(2);

            stress_bot_local(0) = stress_bot_plane(0);
            stress_bot_local(1) = stress_bot_plane(1);
            stress_bot_local(3) = stress_shear(0);
            stress_bot_local(4) = stress_shear(1);
            stress_bot_local(5) = stress_bot_plane(2);

            const Mat3 section_basis = this->section_basis_global(r, s, axes);
            Stress stress_top_global = Stress{stress_top_local}.transform(section_basis);
            Stress stress_bot_global = Stress{stress_bot_local}.transform(section_basis);

            const Index node_idx = static_cast<Index>(node_id);
            for (int j = 0; j < 6; ++j) {
                stress_top(node_idx, j) += stress_top_global(j);
                stress_bot(node_idx, j) += stress_bot_global(j);
            }
        }
    }

    bool supports_shell_resultants() const override { return true; }

    void compute_shell_resultants_nodal(Field& displacement,
                                        Field& resultants) override {
        Mat3 axes = get_xyz_axes();
        Precision h = this->get_section()->thickness_;

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        auto abd = this->get_section()->get_abd();
        Mat2 mat_shear    = this->get_section()->get_shear();

        abd *= topo_scale;
        mat_shear    *= topo_scale;

        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending.setZero();
        disp_shear.setZero();

        for (int a = 0; a < N; ++a) {
            ID node_id = this->nodes()[a];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            disp_membrane(2 * a    ) = disp_xyz(0);
            disp_membrane(2 * a + 1) = disp_xyz(1);

            disp_bending(3 * a    ) = disp_xyz(2);
            disp_bending(3 * a + 1) = disp_rot(0);
            disp_bending(3 * a + 2) = disp_rot(1);

            disp_shear(3 * a    ) = disp_xyz(2);
            disp_shear(3 * a + 1) = disp_rot(0);
            disp_shear(3 * a + 2) = disp_rot(1);
        }

        StaticMatrix<N, 2> coords    = this->geometry.node_coords_local();
        LocalCoords        xy_coords = get_xy_coords(axes);

        for (int i = 0; i < N; ++i) {
            ID node_id = this->nodes()[i];

            Precision r = coords(i, 0);
            Precision s = coords(i, 1);

            ShapeFunction   shape_func = this->shape_function(r, s);
            ShapeDerivative shape_der  = this->shape_derivative(r, s);
            Jacobian        jac        = this->jacobian(shape_der, xy_coords);

            auto B_membrane = this->strain_disp_membrane(shape_der, jac);
            auto B_bending  = this->strain_disp_bending(shape_der, jac);
            auto B_shear    = this->strain_disp_shear_at(r, s, xy_coords);

            const Mat2 A = this->section_to_element_2d(r, s, axes);
            const auto T_plane = this->plane_strain_transform(A);
            const auto T_shear = this->shear_strain_transform(A);

            Vec3 eps_membrane = T_plane * (B_membrane * disp_membrane);
            Vec3 kappa        = T_plane * (B_bending  * disp_bending);
            Vec2 gamma_shear  = T_shear * (B_shear    * disp_shear);

            Vec6 generalized_strain;
            generalized_strain << eps_membrane, kappa;
            Vec6 generalized_resultant = abd * generalized_strain;
            Vec3 membrane_force = generalized_resultant.template segment<3>(0);
            Vec3 bending_moment = generalized_resultant.template segment<3>(3);
            Vec2 shear_force    = mat_shear * gamma_shear;

            const Index node_idx = static_cast<Index>(node_id);
            resultants(node_idx, 0) += membrane_force(0);
            resultants(node_idx, 1) += membrane_force(1);
            resultants(node_idx, 2) += membrane_force(2);
            resultants(node_idx, 3) += bending_moment(0);
            resultants(node_idx, 4) += bending_moment(1);
            resultants(node_idx, 5) += bending_moment(2);
            resultants(node_idx, 6) += shear_force(0);
            resultants(node_idx, 7) += shear_force(1);
        }
    }
    // in DefaultShellElement<...>
    void compute_stress_strain(Field& ip_stress,
                               Field& /*ip_strain*/,
                               Field& displacement,
                               int       ip_offset) override {
        Mat3            axes      = get_xyz_axes();
        auto            xy_coords = get_xy_coords(axes);
        const Precision h         = this->get_section()->thickness_;

        // Membrane material (σx, σy, τxy)
        Mat3 Dm = this->get_section()->get_abd().template block<3, 3>(0, 0) / h;

        // Local membrane displacement vector [ux, uy] per node
        StaticVector<2 * N> u_mem;
        u_mem.setZero();
        for (Index i = 0; i < N; ++i) {
            ID   node_id     = this->nodes()[i];
            Vec6 u_glob      = displacement.row_vec6(static_cast<Index>(node_id));    // (ux,uy,uz,rx,ry,rz)
            Vec3 t_glob      = u_glob.head<3>();
            Vec3 t_loc       = axes * t_glob;
            u_mem(2 * i)     = t_loc(0);    // ux
            u_mem(2 * i + 1) = t_loc(1);    // uy
        }

        const auto& scheme = this->integration_scheme();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const Precision r     = scheme.get_point(ip).r;
            const Precision s     = scheme.get_point(ip).s;

            ShapeDerivative dH_rs = this->shape_derivative(r, s);
            Jacobian        J     = this->jacobian(dH_rs, xy_coords);

            // B-matrix for membrane part (3×2N)
            StaticMatrix<3, 2 * N> Bm = this->strain_disp_membrane(dH_rs, J);

            // Membrane stresses in local coords
            StaticVector<3> sigma = Dm * (Bm * u_mem);

            // Stress resultants (force per unit length)
            const Precision nxx = sigma(0) * h;
            const Precision nyy = sigma(1) * h;
            const Precision nxy = sigma(2) * h;

            const Index     row = static_cast<Index>(ip_offset) + ip;
            ip_stress(row, 0)   = nxx;    // xx
            ip_stress(row, 1)   = nyy;    // yy
            ip_stress(row, 2)   = 0.0;    // zz
            ip_stress(row, 3)   = 0.0;    // yz
            ip_stress(row, 4)   = 0.0;    // zx
            ip_stress(row, 5)   = nxy;    // xy
        }
    }
};
}    // namespace fem::model

#endif    // SHELL_SIMPLE_H
