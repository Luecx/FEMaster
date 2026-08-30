/**
 * @file shell_thermal.h
 * @brief Adds in-plane heat conduction to simple shell elements.
 */

#pragma once

#include "s3.h"
#include "s4.h"
#include "../element/element_thermal.h"

#include <cmath>

namespace fem::model {

template<class Base, Index N>
struct ThermalShell : Base, ThermalElement {
    ThermalShell(ID elem_id, std::array<ID, N> node_ids)
        : Base(elem_id, node_ids) {}

    ElementPtr copy() const override {
        return std::make_shared<ThermalShell<Base, N>>(this->elem_id, this->node_ids);
    }

    MapMatrix conductivity(Precision* buffer) override {
        auto material = this->get_material();
        logging::error(material->has_thermal_conductivity(),
            this->type_name(), ": material has no thermal conductivity for element ", this->elem_id);

        const Precision conductivity = material->get_thermal_conductivity();
        const Precision thickness    = this->get_section()->thickness_;
        auto axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();

        StaticMatrix<N, N> matrix = StaticMatrix<N, N>::Zero();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            auto dH = this->shape_derivative(point.r, point.s);
            auto jac = this->jacobian(dH, xy_coords);
            const Precision det = std::abs(jac.determinant());
            logging::error(det > Precision(0),
                this->type_name(), ": degenerate thermal shell element ", this->elem_id);

            const auto dH_xy = dH * jac.inverse();
            matrix.noalias() += conductivity * thickness * point.w * det
                              * (dH_xy * dH_xy.transpose());
        }

        matrix = Precision(0.5) * (matrix + matrix.transpose());
        MapMatrix mapped(buffer, N, N);
        mapped = matrix;
        return mapped;
    }

    MapMatrix capacity(Precision* buffer) override {
        auto material = this->get_material();
        logging::error(material->has_density(),
            this->type_name(), ": material has no density for element ", this->elem_id);
        logging::error(material->has_thermal_specific_heat(),
            this->type_name(), ": material has no specific heat for element ", this->elem_id);

        const Precision scale = material->get_density()
                              * material->get_thermal_specific_heat()
                              * this->get_section()->thickness_;
        auto axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();

        StaticMatrix<N, N> matrix = StaticMatrix<N, N>::Zero();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            const auto H = this->shape_function(point.r, point.s);
            auto dH = this->shape_derivative(point.r, point.s);
            auto jac = this->jacobian(dH, xy_coords);
            const Precision det = std::abs(jac.determinant());
            matrix.noalias() += scale * point.w * det * (H * H.transpose());
        }

        matrix = Precision(0.5) * (matrix + matrix.transpose());
        MapMatrix mapped(buffer, N, N);
        mapped = matrix;
        return mapped;
    }

    void compute_heat_flux(Field& heat_flux, const Field& temperature) override {
        logging::error(heat_flux.domain == FieldDomain::ELEMENT_IP,
            this->type_name(), ": heat flux output must use ELEMENT_IP domain");
        logging::error(heat_flux.components >= 3,
            this->type_name(), ": heat flux output requires three components");
        logging::error(temperature.domain == FieldDomain::NODE && temperature.components == 1,
            this->type_name(), ": temperature must be a scalar NODE field");

        auto material = this->get_material();
        logging::error(material->has_thermal_conductivity(),
            this->type_name(), ": material has no thermal conductivity for element ", this->elem_id);
        const Precision conductivity = material->get_thermal_conductivity();

        StaticVector<N> local_temperature;
        for (Index i = 0; i < N; ++i) {
            local_temperature(i) = temperature(static_cast<Index>(this->node_ids[i]), 0);
        }

        auto axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);
        const auto& scheme = this->integration_scheme();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            auto dH = this->shape_derivative(point.r, point.s);
            auto jac = this->jacobian(dH, xy_coords);
            const Precision det = std::abs(jac.determinant());
            logging::error(det > Precision(0),
                this->type_name(), ": degenerate thermal shell element ", this->elem_id);

            const auto dH_xy = dH * jac.inverse();
            const Vec2 gradient = dH_xy.transpose() * local_temperature;
            const Vec2 local_flux = -conductivity * gradient;
            const Vec3 global_flux = axes.row(0).transpose() * local_flux(0)
                                   + axes.row(1).transpose() * local_flux(1);

            const Index row = this->ip_index(ip);
            for (Dim component = 0; component < 3; ++component) {
                heat_flux(row, component) = global_flux(component);
            }
        }
    }
};

using ThermalS3 = ThermalShell<S3, 3>;
using ThermalS4 = ThermalShell<S4, 4>;

} // namespace fem::model
