//
// Created by f_eggers on 11.12.2024.
//

#ifndef MITC4_H
#define MITC4_H

#include "shell.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

struct MITC4 : ShellElement<4> {

    Surface4 geometry;
    quadrature::Quadrature integration_scheme_;

    MITC4(ID p_elem_id, const std::array<ID, 4>& p_node_ids)
        : ShellElement<4>(p_elem_id, p_node_ids),
            geometry(p_node_ids),
            integration_scheme_(quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_CUBIC) {
    }
    ~MITC4() override = default;

    Vec4 thicknesses() {
        // TODO
        return {
            this->get_section()->thickness,
            this->get_section()->thickness,
            this->get_section()->thickness,
            this->get_section()->thickness};
    }
    Vec3 normal(Precision r, Precision s) {
        auto node_coords_glob = this->_model_data->get(POSITION);
        return geometry.normal(node_coords_glob, Vec2{r,s});
    }
    StaticMatrix<4, 3> shell_directions() {
        // TODO
        auto               node_coords_local = geometry.node_coords_local();

        Vec3               n1                = normal(node_coords_local(0, 0), node_coords_local(0, 1));
        Vec3               n2                = normal(node_coords_local(1, 0), node_coords_local(1, 1));
        Vec3               n3                = normal(node_coords_local(2, 0), node_coords_local(2, 1));
        Vec3               n4                = normal(node_coords_local(3, 0), node_coords_local(3, 1));

        StaticMatrix<4, 3> res;
        res.row(0) = n1.transpose();
        res.row(1) = n2.transpose();
        res.row(2) = n3.transpose();
        res.row(3) = n4.transpose();
        return res;
    }
    StaticMatrix<4,3> thickness_times_directions() {
        auto shell_dir = shell_directions();
        Vec4 thickness = thicknesses();
        StaticMatrix<4,3> res;
        for(int i = 0; i < 4; i++) {
            res.row(i) = thickness(i) * shell_dir.row(i);
        }
        return res;
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    Mat3 covariant(const Vec3& pos) {
        // derivative of X w.r.t r
        auto node_coords_glob = this->node_coords_global();
        auto dX0drs           = geometry.jacobian(node_coords_glob, pos(0), pos(1));
        auto N                = geometry.shape_function(pos(0), pos(1));
        auto dNdrs            = geometry.shape_derivative(pos(0), pos(1));

        // thickness x normal at each corner (4 x 3)
        auto hv_corner = thickness_times_directions();

        auto dNdr      = dNdrs.col(0);
        auto dNds      = dNdrs.col(1);

        // derivative of (hv(r,s) w.r.t r and s
        auto dNdrhv = hv_corner.transpose() * dNdr;
        auto dNdshv = hv_corner.transpose() * dNds;

        Vec3 dXdr = dX0drs.col(0) + pos(2) / 2 * dNdrhv;
        Vec3 dXds = dX0drs.col(1) + pos(2) / 2 * dNdshv;
        Vec3 dXdt = 0.5 * hv_corner.transpose() * N;

        Mat3 res;
        res << dXdr, dXds, dXdt;
        return res;
    }
    Mat3 contravariant(const Mat3& covariant) {
        Vec3 gr = covariant.col(0);
        Vec3 gs = covariant.col(1);
        Vec3 gt = covariant.col(2);

        Precision grgsgt = gr.dot(gs.cross(gt));

        Vec3 contra_gr = gs.cross(gt) / grgsgt;
        Vec3 contra_gs = gt.cross(gr) / grgsgt;
        Vec3 contra_gt = gr.cross(gs) / grgsgt;

        Mat3 res;
        res << contra_gr, contra_gs, contra_gt;
        return res;
    }
    Mat3 localbasis(const Mat3& covariant) {
        Vec3 gr = covariant.col(0);
        Vec3 gs = covariant.col(1);
        Vec3 gt = covariant.col(2);

        Vec3 e3 = gt / gt.norm();
        Vec3 e1 = gs.cross(e3) / gs.cross(e3).norm();
        Vec3 e2 = e3.cross(e1);
    }

    StaticMatrix<6, 6 * 4> strain_displacement(const Vec3& pos) {
        Precision r = pos(0);
        Precision s = pos(1);
        Precision t = pos(2);

        StaticMatrix<6, 6 * 4> res{};

        // derivatives of shape functions w.r.t r and s
        auto dN = geometry.shape_derivative(r, s);
        auto dNr = dN.col(0);
        auto dNs = dN.col(1);

        auto covariant = this->covariant(pos);
        auto gr = covariant.col(0);
        auto gs = covariant.col(1);
        auto gt = covariant.col(2);

        auto hv_mat = this->thickness_times_directions();

        for(ID local_id = 0; local_id < 4; local_id++) {

            Index dof_1 = 6 * local_id;
            Index dof_2 = 6 * local_id + 1;
            Index dof_3 = 6 * local_id + 2;
            Index dof_4 = 6 * local_id + 3;
            Index dof_5 = 6 * local_id + 4;
            Index dof_6 = 6 * local_id + 5;

            // extract variables for this dof
            Vec3 hv = hv_mat.col(local_id);

            auto dNi_dr = dNr(local_id);
            auto dNi_ds = dNs(local_id);

            int x = 0;
            int y = 1;
            int z = 2;

            // strain e11
            // u_1x  gx*Derivative(N1(r, s), r)
            // u_1y  gy*Derivative(N1(r, s), r)
            // u_1z  gz*Derivative(N1(r, s), r)
            // theta_1x -gy*h_1*t*v_1z*Derivative(N1(r, s), r)/2 + gz*h_1*t*v_1y*Derivative(N1(r, s), r)/2
            // theta_1y  gx*h_1*t*v_1z*Derivative(N1(r, s), r)/2 - gz*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            // theta_1z -gx*h_1*t*v_1y*Derivative(N1(r, s), r)/2 + gy*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            res(0, dof_1) = gr(x) * dNi_dr;
            res(0, dof_2) = gr(y) * dNi_dr;
            res(0, dof_3) = gr(z) * dNi_dr;

            res(0, dof_4) = (gr(z) * hv(y) - gr(y) * hv(z)) * t * dNi_dr / 2;
            res(0, dof_5) = (gr(x) * hv(z) - gr(z) * hv(x)) * t * dNi_dr / 2;
            res(0, dof_6) = (gr(y) * hv(x) - gr(x) * hv(y)) * t * dNi_dr / 2;

            // strain e22
            // u_1x gx*Derivative(N1(r, s), s)
            // u_1y gy*Derivative(N1(r, s), s)
            // u_1z gz*Derivative(N1(r, s), s)
            // theta_1x -gy*h_1*t*v_1z*Derivative(N1(r, s), s)/2 + gz*h_1*t*v_1y*Derivative(N1(r, s), s)/2
            // theta_1y  gx*h_1*t*v_1z*Derivative(N1(r, s), s)/2 - gz*h_1*t*v_1x*Derivative(N1(r, s), s)/2
            // theta_1z -gx*h_1*t*v_1y*Derivative(N1(r, s), s)/2 + gy*h_1*t*v_1x*Derivative(N1(r, s), s)/2
            res(1, dof_1) = gs(x) * dNi_ds;
            res(1, dof_2) = gs(y) * dNi_ds;
            res(1, dof_3) = gs(z) * dNi_ds;
            res(1, dof_4) = (gs(z) * hv(y) - gs(y) * hv(z)) * t * dNi_ds / 2;
            res(1, dof_5) = (gs(x) * hv(z) - gs(z) * hv(x)) * t * dNi_ds / 2;
            res(1, dof_6) = (gs(y) * hv(x) - gs(x) * hv(y)) * t * dNi_ds / 2;

            // strain e33 = 0
            res(2, dof_1) = 0;
            res(2, dof_2) = 0;
            res(2, dof_3) = 0;
            res(2, dof_4) = 0;
            res(2, dof_5) = 0;
            res(2, dof_6) = 0;

            // strain e23
            // u_1x gtx*Derivative(N1(r, s), s)
            // u_1y gty*Derivative(N1(r, s), s)
            // u_1z gtz*Derivative(N1(r, s), s)
            // theta_1x -gsy*h_1*v_1z*N1(r, s)/2 + gsz*h_1*v_1y*N1(r, s)/2 - gty*h_1*t*v_1z*Derivative(N1(r, s), s)/2 + gtz*h_1*t*v_1y*Derivative(N1(r, s), s)/2
            // theta_1y gsx*h_1*v_1z*N1(r, s)/2 - gsz*h_1*v_1x*N1(r, s)/2 + gtx*h_1*t*v_1z*Derivative(N1(r, s), s)/2 - gtz*h_1*t*v_1x*Derivative(N1(r, s), s)/2
            // theta_1z -gsx*h_1*v_1y*N1(r, s)/2 + gsy*h_1*v_1x*N1(r, s)/2 - gtx*h_1*t*v_1y*Derivative(N1(r, s), s)/2 + gty*h_1*t*v_1x*Derivative(N1(r, s), s)/2
            res(3, dof_1) = gt(x) * dNi_ds;
            res(3, dof_2) = gt(y) * dNi_ds;
            res(3, dof_3) = gt(z) * dNi_ds;
            res(3, dof_4) = (gs(z) * hv(y) - gs(y) * hv(z)) * dNi_ds / 2;
            res(3, dof_5) = (gs(x) * hv(z) - gs(z) * hv(x)) * dNi_ds / 2;
            res(3, dof_6) = (gs(y) * hv(x) - gs(x) * hv(y)) * dNi_ds / 2;

            //---- e13 ----
            // u_1x gtx*Derivative(N1(r, s), r)
            // u_1y gty*Derivative(N1(r, s), r)
            // u_1z gtz*Derivative(N1(r, s), r)
            // theta_1x -gry*h_1*v_1z*N1(r, s)/2 + grz*h_1*v_1y*N1(r, s)/2 - gty*h_1*t*v_1z*Derivative(N1(r, s), r)/2 + gtz*h_1*t*v_1y*Derivative(N1(r, s), r)/2
            // theta_1y grx*h_1*v_1z*N1(r, s)/2 - grz*h_1*v_1x*N1(r, s)/2 + gtx*h_1*t*v_1z*Derivative(N1(r, s), r)/2 - gtz*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            // theta_1z -grx*h_1*v_1y*N1(r, s)/2 + gry*h_1*v_1x*N1(r, s)/2 - gtx*h_1*t*v_1y*Derivative(N1(r, s), r)/2 + gty*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            res(4, dof_1) = gt(x) * dNi_dr;
            res(4, dof_2) = gt(y) * dNi_dr;
            res(4, dof_3) = gt(z) * dNi_dr;
            res(4, dof_4) = (gr(z) * hv(y) - gr(y) * hv(z)) * dNi_dr / 2;
            res(4, dof_5) = (gr(x) * hv(z) - gr(z) * hv(x)) * dNi_dr / 2;
            res(4, dof_6) = (gr(y) * hv(x) - gr(x) * hv(y)) * dNi_dr / 2;

            // u_1x grx*Derivative(N1(r, s), s) + gsx*Derivative(N1(r, s), r)
            // u_1y gry*Derivative(N1(r, s), s) + gsy*Derivative(N1(r, s), r)
            // u_1z grz*Derivative(N1(r, s), s) + gsz*Derivative(N1(r, s), r)
            // theta_1x -gry*h_1*t*v_1z*Derivative(N1(r, s), s)/2 + grz*h_1*t*v_1y*Derivative(N1(r, s), s)/2 - gsy*h_1*t*v_1z*Derivative(N1(r, s), r)/2 + gsz*h_1*t*v_1y*Derivative(N1(r, s), r)/2
            // theta_1y grx*h_1*t*v_1z*Derivative(N1(r, s), s)/2 - grz*h_1*t*v_1x*Derivative(N1(r, s), s)/2 + gsx*h_1*t*v_1z*Derivative(N1(r, s), r)/2 - gsz*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            // theta_1z -grx*h_1*t*v_1y*Derivative(N1(r, s), s)/2 + gry*h_1*t*v_1x*Derivative(N1(r, s), s)/2 - gsx*h_1*t*v_1y*Derivative(N1(r, s), r)/2 + gsy*h_1*t*v_1x*Derivative(N1(r, s), r)/2
            res(5, dof_1) = gr(x) * dNi_ds + gs(x) * dNi_dr;
            res(5, dof_2) = gr(y) * dNi_ds + gs(y) * dNi_dr;
            res(5, dof_3) = gr(z) * dNi_ds + gs(z) * dNi_dr;
            res(5, dof_4) = (gr(z) * hv(y) - gr(y) * hv(z)) * dNi_ds / 2 + (gs(z) * hv(y) - gs(y) * hv(z)) * dNi_dr / 2;
            res(5, dof_5) = (gr(x) * hv(z) - gr(z) * hv(x)) * dNi_ds / 2 + (gs(x) * hv(z) - gs(z) * hv(x)) * dNi_dr / 2;
            res(5, dof_6) = (gr(y) * hv(x) - gr(x) * hv(y)) * dNi_ds / 2 + (gs(y) * hv(x) - gs(x) * hv(y)) * dNi_dr / 2;
        }
    }


    SurfacePtr surface(ID surface_id) override {
        (void) surface_id;
        return nullptr;
    }
    Precision  volume() override {
        return 0;
    }
    MapMatrix  stiffness(Precision* buffer) override {
        auto covariant = this->covariant(Vec3{0,0,0});
        std::cout << covariant << std::endl;
        exit(0);
        return MapMatrix(buffer, 4, 4);
    }
    MapMatrix  mass(Precision* buffer) override {
        return MapMatrix(buffer, 4, 4);
    }
};

}





#endif //MITC4_H
