//
// Created by f_eggers on 11.12.2024.
//

#ifndef S4_H
#define S4_H

#include "shell.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface8.h"
#include "../geometry/surface/surface3.h"
#include "../geometry/surface/surface6.h"

// this file is partially based on the implementation of
// https://github.com/JWock82/Pynite/blob/main/Archived/S4.py
//

namespace fem::model {

template<Index N, typename SFType, quadrature::Domain INT_D, quadrature::Order INT_O>
struct DefaultShellElement : public ShellElement<N> {
    SFType geometry;
    quadrature::Quadrature integration_scheme_;

    DefaultShellElement(ID p_elem_id, std::array<ID, N> p_node)
        : ShellElement<N>(p_elem_id, p_node)
        , geometry(p_node)
        , integration_scheme_(INT_D, INT_O) {}
    ~DefaultShellElement() override = default;

    // when assuming a planar element, the normal is constant
    // we can also project the shell and use a local x,y system
    // this must not be confused with the local r,s system which is used for integration
    Mat3 get_xyz_axes() {
        StaticMatrix<N, 3> node_coords = this->node_coords_global();

        Vec3 n1 = node_coords.row(0);
        Vec3 n2 = node_coords.row(1);
        Vec3 n3 = node_coords.row(2);

        Vec3 n12 = n2 - n1;
        Vec3 n13 = n3 - n1;

        Vec3 x_axis = n12 / n12.norm();
        Vec3 z_axis = n12.cross(n13) / (n12.cross(n13)).norm();
        Vec3 y_axis = z_axis.cross(x_axis);

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
    StaticMatrix<N*6,N*6> transformation() {
        auto axes = get_xyz_axes();
        StaticMatrix<N*6,N*6> res = StaticMatrix<N*6,N*6>::Zero();
        for (int i = 0; i < 2*N; ++i) { // Loop over the 8 nodes
            res.template block<3, 3>(3 * i, 3 * i) = axes;
        }
        return res;
    }


    StaticMatrix<N, 2> get_xy_coords() {
        StaticMatrix<N, 3> node_coords = this->node_coords_global();

        Mat3 axes = get_xyz_axes();
        Vec3 x_axis = axes.row(0).transpose();
        Vec3 y_axis = axes.row(1).transpose();

        // for each node calculate the local x,y coordinates by projecting onto the axes
        StaticMatrix<N, 2> xy_coords;
        for (Index i = 0; i < N; i++) {
            Vec3 ni = node_coords.row(i);
            Vec3 n0 = node_coords.row(0);
            Vec3 n = ni - n0;
            xy_coords(i, 0) = n.dot(x_axis);
            xy_coords(i, 1) = n.dot(y_axis);
        }
        return xy_coords;
    }

    // when assuming a planar element, the normal is constant
    // one can simply retrieve it from the local xyz axes
    Vec3 get_normal() {
        return get_xyz_axes().row(2).transpose();
    }

    auto shape_function(Precision r, Precision s) {
        return geometry.shape_function(r, s);
    }

    auto shape_derivative(Precision r, Precision s) {
        return geometry.shape_derivative(r, s);
    }

    // function which relates derivatives in the local r,s to the local x,y system
    // the entries will contain the mapping
    // [dx/dr, dy/dr]
    // [dx/ds, dy/ds]
    // x = N1*x1 + N2*x2 + N3*x3 + N4*x4
    // y = N1*y1 + N2*y2 + N3*y3 + N4*y4
    // hence the derivatives are
    // dx/dr = dN1/dr*x1 + dN2/dr*x2 + dN3/dr*x3 + dN4/dr*x4
    Mat2 jacobian() {
        // computing the jacobian is done by computing the derivatives of the shape functions
        // with the local x,y system
        auto xy_coords = get_xy_coords();
        auto shape_derivatives = shape_derivative(0, 0);

        Mat2 jacobian;
        jacobian(0, 0) = shape_derivatives.col(0).dot(xy_coords.col(0));
        jacobian(0, 1) = shape_derivatives.col(0).dot(xy_coords.col(1));
        jacobian(1, 0) = shape_derivatives.col(1).dot(xy_coords.col(0));
        jacobian(1, 1) = shape_derivatives.col(1).dot(xy_coords.col(1));
        return jacobian;
    }


    StaticMatrix<3, N*3> strain_disp_bending(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();

        Mat2 inv = jacobian.inverse();

        auto dH = (inv * shape_der.transpose());


        StaticMatrix<3, N*3> res{};

        // dofs are displacement in z, rotation around x, rotation around y
        // its only a function of the rotational dofs (second derivative of displacement)
        for(int i = 0; i < N; i++) {
            res(0, 3*i  ) = 0;
            res(0, 3*i+1) = 0;
            res(0, 3*i+2) = -dH(0, i);

            res(1, 3*i  ) = 0;
            res(1, 3*i+1) = dH(1, i);
            res(1, 3*i+2) = 0;

            res(2, 3*i  ) = 0;
            res(2, 3*i+1) = dH(0, i);
            res(2, 3*i+2) = -dH(1, i);
        }
        // res <<
        //     0,    0,     -dH(0, 0), 0,    0,     -dH(0, 1), 0,    0,     -dH(0, 2), 0,    0,     -dH(0, 3),
        //     0, dH(1, 0),     0,     0, dH(1, 1),     0,     0, dH(1, 2),     0,     0, dH(1, 3),     0    ,
        //     0, dH(0, 0), -dH(1, 0), 0, dH(0, 1), -dH(1, 1), 0, dH(0, 2), -dH(1, 2), 0, dH(0, 3), -dH(1, 3);
        return res;
    }

    StaticMatrix<2, 3 * N> strain_disp_shear(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();
        auto H = shape_function(r, s);
        auto dH = (jacobian.inverse() * shape_der.transpose());

        // StaticMatrix<2, 12> res{};
        // // dofs are displacement in z, rotation around x, rotation around y
        // // its only a function of the rotational dofs (second derivative of displacement)
        // res <<
        //     dH(0,0),  0, H(0), dH(0,1),  0, H(1), dH(0,2),  0, H(2), dH(0,3),  0, H(3),
        //     dH(1,0),-H(0),  0, dH(1,1),-H(1),  0, dH(1,2),-H(2),  0, dH(1,3),-H(3), 0;
        // return res;
        StaticMatrix<2, 3 * N> res{};
        for(int i = 0; i < N; i++) {
            res(0, 3*i  ) = dH(0, i);
            res(0, 3*i+1) = 0;
            res(0, 3*i+2) = H(i);

            res(1, 3*i  ) = dH(1, i);
            res(1, 3*i+1) = -H(i);
            res(1, 3*i+2) = 0;
        }
        return res;
    }

    StaticMatrix<3, 2 * N> strain_disp_membrane(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();

        auto dH = (jacobian.inverse() * shape_der.transpose());
        // StaticMatrix<3, 8> res{};
        // res << dH(0,0), 0      , dH(0,1), 0      , dH(0,2), 0      , dH(0,3), 0      ,
        //            0  , dH(1,0), 0      , dH(1,1), 0      , dH(1,2), 0      , dH(1,3),
        //        dH(1,0), dH(0,0), dH(1,1), dH(0,1), dH(1,2), dH(0,2), dH(1,3), dH(0,3);
        StaticMatrix<3, 2 * N> res{};
        for(int i = 0; i < N; i++) {
            res(0, 2*i  ) = dH(0, i);
            res(0, 2*i+1) = 0;

            res(1, 2*i  ) = 0;
            res(1, 2*i+1) = dH(1, i);

            res(2, 2*i  ) = dH(1, i);
            res(2, 2*i+1) = dH(0, i);
        }
        return res;
    }

    StaticMatrix<6 * N, 6 * N> stiffness_bending() {

        auto mat_bend  = this->get_elasticity()->get_bend(this->get_section()->thickness);
        auto mat_shear = this->get_elasticity()->get_shear(this->get_section()->thickness);

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_bend =
            [this, mat_bend](Precision r, Precision s, Precision t) -> StaticMatrix<3 * N, 3 * N> {
                Precision det;

                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_bending(r, s);
                auto E = mat_bend;
                return B.transpose() * (E * B) * jac_det;
        };

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_shear =
            [this, mat_shear](Precision r, Precision s, Precision t) -> StaticMatrix<3 * N, 3 * N> {
                Precision det;
                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_shear(r, s);
                auto E = mat_shear;
                return B.transpose() * (E * B) * jac_det;
        };


        StaticMatrix<3 * N, 3 * N> stiff_bend  = this->integration_scheme().integrate(func_bend);
        StaticMatrix<3 * N, 3 * N> stiff_shear = this->integration_scheme().integrate(func_shear);
        StaticMatrix<3 * N, 3 * N> stiff_local = stiff_bend + stiff_shear;

        Precision k_drill = 1e32;
        for(int i = 0; i < N; i++) {
            k_drill = std::min(k_drill, stiff_local(3*i+1, 3*i+1));
            k_drill = std::min(k_drill, stiff_local(3*i+2, 3*i+2));
        }
        k_drill /= 1000;

        StaticMatrix<6*N, 6*N> res;
        res.setZero();

        int dof_map[] {
            2, // z
            3, // rx
            4, // ry
        };

        for(Index i = 0; i < 3*N; i++) {
            for(Index j = 0; j < 3*N; j++) {
                auto i_local_id = i / 3;
                auto j_local_id = j / 3;
                auto i_dof     = i % 3;
                auto j_dof     = j % 3;

                auto i_glob_index = i_local_id * 6 + dof_map[i_dof];
                auto j_glob_index = j_local_id * 6 + dof_map[j_dof];

                res(i_glob_index, j_glob_index) = stiff_local(i, j);
            }
        }

        for(int i = 0; i < N; i++) {
            res(6*i+5, 6*i+5) = k_drill;
        }

        return res;
    }

    StaticMatrix<6*N, 6*N> stiffness_membrane() {
        auto mat_membrane = this->get_elasticity()->get_memb();
        std::function<StaticMatrix<2*N, 2*N>(Precision, Precision, Precision)> func_membrane =
            [this, &mat_membrane](Precision r, Precision s, Precision t) -> StaticMatrix<2*N, 2*N> {
                Precision det;
                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_membrane(r, s);
                auto E = mat_membrane;
                return B.transpose() * (E * B) * jac_det * this->get_section()->thickness;
        };

        StaticMatrix<2*N, 2*N> stiff_membrane = this->integration_scheme().integrate(func_membrane);
        StaticMatrix<6*N, 6*N> res;
        res.setZero();

        int dof_map [] {
            0, 1
        };
        for(Index i = 0; i < 2*N; i++) {
            for(Index j = 0; j < 2*N; j++) {
                auto i_local_id = i / 2;
                auto j_local_id = j / 2;
                auto i_dof     = i % 2;
                auto j_dof     = j % 2;

                auto i_glob_index = i_local_id * 6 + dof_map[i_dof];
                auto j_glob_index = j_local_id * 6 + dof_map[j_dof];

                res(i_glob_index, j_glob_index) = stiff_membrane(i, j);
            }
        }
        return res;
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }


    SurfacePtr surface(ID surface_id) override {
        (void) surface_id;
        return nullptr;
    }
    Precision  volume() override {
        return 0;
    }
    MapMatrix  stiffness(Precision* buffer) override {
        MapMatrix mapped{buffer, 6*N, 6*N};
        auto trans = transformation();
        auto stiff = stiffness_bending() + stiffness_membrane();
        mapped = trans.transpose() * stiff * trans;
        return mapped;
    }
    MapMatrix  mass(Precision* buffer) override {
        return MapMatrix(buffer, 4, 4);
    }


};


struct S4 : DefaultShellElement<4, Surface4, quadrature::Domain::DOMAIN_ISO_QUAD, quadrature::Order::ORDER_CUBIC> {
    S4(ID p_elem_id, std::array<ID, 4> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}
};
struct S8 : DefaultShellElement<8, Surface8, quadrature::Domain::DOMAIN_ISO_QUAD, quadrature::Order::ORDER_QUINTIC> {
    S8(ID p_elem_id, std::array<ID, 8> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}
};
struct S3 : DefaultShellElement<3, Surface3, quadrature::Domain::DOMAIN_ISO_TRI, quadrature::Order::ORDER_QUINTIC> {
    S3(ID p_elem_id, std::array<ID, 3> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}
};
struct S6 : DefaultShellElement<6, Surface6, quadrature::Domain::DOMAIN_ISO_TRI, quadrature::Order::ORDER_QUINTIC> {
    S6(ID p_elem_id, std::array<ID, 6> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}
};


}





#endif //S4_H
