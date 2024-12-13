//
// Created by f_eggers on 11.12.2024.
//

#ifndef S4_H
#define S4_H

#include "shell.h"
#include "../geometry/surface/surface4.h"

// this file is partially based on the implementation of
// https://github.com/JWock82/Pynite/blob/main/Archived/S4.py
//

namespace fem::model {

struct S4 : ShellElement<4> {

    Surface4 geometry;
    quadrature::Quadrature integration_scheme_;

    S4(ID p_elem_id, const std::array<ID, 4>& p_node_ids)
        : ShellElement<4>(p_elem_id, p_node_ids),
            geometry(p_node_ids),
            integration_scheme_(quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_CUBIC) {
    }
    ~S4() override = default;


    // when assuming a planar element, the normal is constant
    // we can also project the shell and use a local x,y system
    // this must not be confused with the local r,s system which is used for integration
    Mat3 get_xyz_axes() {
        StaticMatrix<4, 3> node_coords = this->node_coords_global();
        StaticMatrix<4, 2> xy_coords;

        Vec3 n1 = node_coords.row(0);
        Vec3 n2 = node_coords.row(1);
        Vec3 n3 = node_coords.row(2);
        Vec3 n4 = node_coords.row(3);

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

    StaticMatrix<24,24> transformation() {
        // Assumption is that the element is planar
        // deviations of that will produce inaccurate results
        // finding the transformation matrix is straight forward
        // the x-axis will be aligned with the direction from node 1 to node 2
        // the z-axis will be perpendicular to the plane of the element
        // the y-axis will be the cross product of the other two
        auto axes = get_xyz_axes();
        StaticMatrix<24,24> res = StaticMatrix<24,24>::Zero();
        for (int i = 0; i < 8; ++i) { // Loop over the 8 nodes
            res.block<3, 3>(3 * i, 3 * i) = axes;
        }
        return res;
    }

    // when assuming a planar element, the normal is constant
    // we can also project the shell and use a local x,y system
    // this must not be confused with the local r,s system which is used for integration
    Mat42 get_xy_coords() {
        StaticMatrix<4, 3> node_coords = this->node_coords_global();
        StaticMatrix<4, 2> xy_coords;

        Vec3 n1 = node_coords.row(0);
        Vec3 n2 = node_coords.row(1);
        Vec3 n3 = node_coords.row(2);
        Vec3 n4 = node_coords.row(3);

        Vec3 n12 = n2 - n1;
        Vec3 n13 = n3 - n1;
        Vec3 n14 = n4 - n1;

        Vec3 x_axis = n12 / n12.norm();
        Vec3 z_axis = n12.cross(n13) / (n12.cross(n13)).norm();
        Vec3 y_axis = z_axis.cross(x_axis);

        // for each node calculate the local x,y coordinates by projecting onto the axes
        xy_coords(0, 0) = 0;
        xy_coords(0, 1) = 0;

        xy_coords(1, 0) = n12.dot(x_axis);
        xy_coords(1, 1) = n12.dot(y_axis);

        xy_coords(2, 0) = n13.dot(x_axis);
        xy_coords(2, 1) = n13.dot(y_axis);

        xy_coords(3, 0) = n14.dot(x_axis);
        xy_coords(3, 1) = n14.dot(y_axis);
        return xy_coords;
    }

    // when assuming a planar element, the normal is constant
    // one can simply retrieve it from the local xyz axes
    Vec3 get_normal() {
        return get_xyz_axes().row(2).transpose();
    }

    Vec4 shape_function(Precision r, Precision s) {
        return geometry.shape_function(r, s);
    }

    Mat42 shape_derivative(Precision r, Precision s) {
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
        Mat42 xy_coords = get_xy_coords();
        Mat42 shape_derivatives = shape_derivative(0, 0);

        Mat2 jacobian;
        jacobian(0, 0) = shape_derivatives.col(0).dot(xy_coords.col(0));
        jacobian(0, 1) = shape_derivatives.col(0).dot(xy_coords.col(1));
        jacobian(1, 0) = shape_derivatives.col(1).dot(xy_coords.col(0));
        jacobian(1, 1) = shape_derivatives.col(1).dot(xy_coords.col(1));
        return jacobian;
    }

    StaticMatrix<3, 12> strain_disp_bending(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();

        Mat2 inv = jacobian.inverse();

        auto dH = (inv * shape_der.transpose());


        StaticMatrix<3, 12> res{};
        // dofs are displacement in z, rotation around x, rotation around y
        // its only a function of the rotational dofs (second derivative of displacement)
        res <<
            0,    0,     -dH(0, 0), 0,    0,     -dH(0, 1), 0,    0,     -dH(0, 2), 0,    0,     -dH(0, 3),
            0, dH(1, 0),     0,     0, dH(1, 1),     0,     0, dH(1, 2),     0,     0, dH(1, 3),     0    ,
            0, dH(0, 0), -dH(1, 0), 0, dH(0, 1), -dH(1, 1), 0, dH(0, 2), -dH(1, 2), 0, dH(0, 3), -dH(1, 3);
        return res;
    }

    StaticMatrix<2, 12> strain_disp_shear(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();
        auto H = shape_function(r, s);
        auto dH = (jacobian.inverse() * shape_der.transpose());

        StaticMatrix<2, 12> res{};
        // dofs are displacement in z, rotation around x, rotation around y
        // its only a function of the rotational dofs (second derivative of displacement)
        res <<
            dH(0,0),  0, H(0), dH(0,1),  0, H(1), dH(0,2),  0, H(2), dH(0,3),  0, H(3),
            dH(1,0),-H(0),  0, dH(1,1),-H(1),  0, dH(1,2),-H(2),  0, dH(1,3),-H(3), 0;
        return res;
    }

    StaticMatrix<3, 8> strain_disp_membrane(Precision r, Precision s) {
        auto shape_der = shape_derivative(r, s);
        auto jacobian = this->jacobian();

        auto dH = (jacobian.inverse() * shape_der.transpose());
        StaticMatrix<3, 8> res{};
        res << dH(0,0), 0      , dH(0,1), 0      , dH(0,2), 0      , dH(0,3), 0      ,
                   0  , dH(1,0), 0      , dH(1,1), 0      , dH(1,2), 0      , dH(1,3),
               dH(1,0), dH(0,0), dH(1,1), dH(0,1), dH(1,2), dH(0,2), dH(1,3), dH(0,3);
        return res;
    }

    StaticMatrix<24, 24> stiffness_bending() {
        auto mat_bend  = get_elasticity()->get_bend(this->get_section()->thickness);
        auto mat_shear = get_elasticity()->get_shear(this->get_section()->thickness);

        std::function<StaticMatrix<12, 12>(Precision, Precision, Precision)> func_bend =
            [this, mat_bend](Precision r, Precision s, Precision t) -> StaticMatrix<12, 12> {
                Precision det;

                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_bending(r, s);
                auto E = mat_bend;
                return B.transpose() * (E * B) * jac_det;
        };

        std::function<StaticMatrix<12, 12>(Precision, Precision, Precision)> func_shear =
            [this, mat_shear](Precision r, Precision s, Precision t) -> StaticMatrix<12, 12> {
                Precision det;
                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_shear(r, s);
                auto E = mat_shear;
                return B.transpose() * (E * B) * jac_det;
        };


        StaticMatrix<12, 12> stiff_bend = this->integration_scheme().integrate(func_bend);
        StaticMatrix<12, 12> stiff_shear = this->integration_scheme().integrate(func_shear);
        StaticMatrix<12, 12> stiff_local = stiff_bend + stiff_shear;

        std::cout << stiff_local << std::endl;

        Precision            k_drill     = std::min({stiff_local(1, 1),
                                                     stiff_local(2, 2),
                                                     stiff_local(4, 4),
                                                     stiff_local(5, 5),
                                                     stiff_local(7, 7),
                                                     stiff_local(8, 8),
                                                     stiff_local(10, 10),
                                                     stiff_local(11, 11)}) / 1000;
        StaticMatrix<24, 24> res;
        res.setZero();

        int index_map[] {
            2,  // z
            3,  // rx
            4,  // ry
            8,  // z
            9,  // rx
            10, // ry
            14, // z
            15, // rx
            16, // ry
            20, // z
            21, // rx
            22  // ry
        };

        for(Index i = 0; i < 12; i++) {
            for(Index j = 0; j < 12; j++) {
                res(index_map[i], index_map[i]) = stiff_local(i, j);
            }
        }
        res(5,5) = k_drill;
        res(11,11) = k_drill;
        res(17,17) = k_drill;
        res(23,23) = k_drill;

        return res;
    }

    StaticMatrix<24, 24> stiffness_membrane() {
        auto mat_membrane = get_elasticity()->get_memb();
        std::function<StaticMatrix<8, 8>(Precision, Precision, Precision)> func_membrane =
            [this, &mat_membrane](Precision r, Precision s, Precision t) -> StaticMatrix<8, 8> {
                Precision det;
                auto jac = this->jacobian();
                auto jac_det = jac.determinant();
                auto B = this->strain_disp_membrane(r, s);
                auto E = mat_membrane;
                return B.transpose() * (E * B) * jac_det;
        };

        int index_map[] {
            0, 1,
            6, 7,
            12, 13,
            18, 19
        };
        StaticMatrix<8, 8> stiff_membrane = this->integration_scheme().integrate(func_membrane);
        std::cout << stiff_membrane << std::endl;

        StaticMatrix<24, 24> res;
        res.setZero();
        for(Index i = 0; i < 8; i++) {
            for(Index j = 0; j < 8; j++) {
                res(index_map[i], index_map[i]) = stiff_membrane(i, j);
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
        MapMatrix mapped{buffer, 24, 24};
        auto trans = transformation();

        std::cout << trans << std::endl;
        try {
            std::cout << stiffness_bending() << std::endl;
        } catch (std::exception& e) {
            std::cout << e.what() << std::endl;
        }
        try {
            std::cout << stiffness_membrane() << std::endl;
        } catch (std::exception& e) {
            std::cout << e.what() << std::endl;
        }

        auto stiff = stiffness_bending() + stiffness_membrane();
        mapped = trans.transpose() * stiff * trans;
        return mapped;
    }
    MapMatrix  mass(Precision* buffer) override {
        return MapMatrix(buffer, 4, 4);
    }




};

}





#endif //S4_H
