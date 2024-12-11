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
        // auto node_coords_glob = this->node_coords_global();
        // auto dX0drs = geometry.jacobian(node_coords_glob, pos(0), pos(1));
        // auto N      = geometry.shape_function(pos(0), pos(1));
        // auto dNdrs  = geometry.shape_derivative(pos(0), pos(1));
        // auto hv     = thickness_times_directions();
        //
        // auto dNdr   = dNdrs.col(0);
        // auto dNds   = dNdrs.col(1);
        //
        // auto dNdrhv = dNdr.asDiagonal() * hv;
        // auto dNdshv = dNds.asDiagonal() * hv;
        //
        // std::cout << dX0drs << std::endl;

        // Vec3 dX0dr = dX0drs.col(0) + pos(2) / 2 *  dNdrhv;
        // Vec3 dX0ds = dX0drs.col(1) + pos(2) / 2 *  dNdshv;
        //
        // Vec3 dXdt  = 0.5 * geometry.interpolate(hv, pos(0), pos(1));
        return {};
        // return geometry.covariant(pos);
    }
    Mat3 contravariant(const Vec3& pos) {
        return {};
            // return geometry.contravariant(pos);
    }

    SurfacePtr surface(ID surface_id) override {
        return nullptr;
    }
    Precision  volume() override {
        return 0;
    }
    MapMatrix  stiffness(Precision* buffer) override {
        auto covariant = this->covariant(Vec3{0,0,0});
        std::cout << covariant << std::endl;
        return MapMatrix(buffer, 4, 4);
        exit(0);
    }
    MapMatrix  mass(Precision* buffer) override {
        return MapMatrix(buffer, 4, 4);
    }
};

}





#endif //MITC4_H
