/**
 * @file b31.h
 * @brief Declares a two-node beam with axial heat-conduction capability.
 */

#pragma once

#include "b33.h"
#include "../element/element_thermal.h"

namespace fem::model {

/**
 * @brief FEMaster B31 compatibility element using the existing B33 mechanics.
 *
 * The thermal formulation assumes temperature is uniform over the beam cross
 * section and conducts only along the reference centerline.
 */
struct B31 : B33, ThermalElement {
    B31(ID elem_id, std::array<ID, 2> node_ids)
        : B33(elem_id, node_ids) {}

    ElementPtr copy() const override {
        return std::make_shared<B31>(elem_id, node_ids);
    }

    std::string type_name() const override { return "B31"; }

    Precision thermal_length() const {
        logging::error(_model_data != nullptr && _model_data->positions_reference != nullptr,
            "B31: reference positions are not initialized for element ", elem_id);
        const Vec3 x0 = _model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 x1 = _model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[1]));
        const Precision length = (x1 - x0).norm();
        logging::error(length > Precision(0),
            "B31: zero reference length for element ", elem_id);
        return length;
    }

    Vec3 thermal_direction() const {
        const Vec3 x0 = _model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 x1 = _model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[1]));
        return (x1 - x0) / thermal_length();
    }

    MapMatrix conductivity(Precision* buffer) override {
        auto material = get_material();
        logging::error(material->has_thermal_conductivity(),
            "B31: material has no thermal conductivity for element ", elem_id);

        const Precision scale = material->get_thermal_conductivity()
                              * get_profile()->area_ / thermal_length();

        MapMatrix matrix(buffer, 2, 2);
        matrix << scale, -scale,
                 -scale,  scale;
        return matrix;
    }

    MapMatrix capacity(Precision* buffer) override {
        auto material = get_material();
        logging::error(material->has_density(),
            "B31: material has no density for element ", elem_id);
        logging::error(material->has_thermal_specific_heat(),
            "B31: material has no specific heat for element ", elem_id);

        const Precision scale = material->get_density()
                              * material->get_thermal_specific_heat()
                              * get_profile()->area_ * thermal_length() / Precision(6);

        MapMatrix matrix(buffer, 2, 2);
        matrix << Precision(2) * scale, scale,
                  scale, Precision(2) * scale;
        return matrix;
    }

    void compute_heat_flux(Field& heat_flux, const Field& temperature) override {
        logging::error(heat_flux.domain == FieldDomain::ELEMENT_IP,
            "B31: heat flux output must use ELEMENT_IP domain");
        logging::error(heat_flux.components >= 3,
            "B31: heat flux output requires three components");
        logging::error(temperature.domain == FieldDomain::NODE && temperature.components == 1,
            "B31: temperature must be a scalar NODE field");

        auto material = get_material();
        logging::error(material->has_thermal_conductivity(),
            "B31: material has no thermal conductivity for element ", elem_id);

        const Precision t0 = temperature(static_cast<Index>(node_ids[0]), 0);
        const Precision t1 = temperature(static_cast<Index>(node_ids[1]), 0);
        const Vec3 flux = -material->get_thermal_conductivity()
                        * (t1 - t0) / thermal_length() * thermal_direction();

        const Index row = ip_index(0);
        for (Dim component = 0; component < 3; ++component) {
            heat_flux(row, component) = flux(component);
        }
    }
};

} // namespace fem::model
