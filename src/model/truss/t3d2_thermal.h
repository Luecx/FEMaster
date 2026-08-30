/**
 * @file t3d2_thermal.h
 * @brief Adds one-dimensional heat conduction to the two-node truss topology.
 */

#pragma once

#include "truss.h"
#include "../element/element_thermal.h"

namespace fem::model {

struct ThermalT3D2 : T3, ThermalElement {
    ThermalT3D2(ID elem_id, std::array<ID, 2> node_ids)
        : T3(elem_id, node_ids) {}

    ElementPtr copy() const override {
        return std::make_shared<ThermalT3D2>(elem_id, node_ids);
    }

    std::string type_name() const override { return "T3D2"; }

    MapMatrix conductivity(Precision* buffer) override {
        auto material = get_material();
        logging::error(material->has_thermal_conductivity(),
            "T3D2: material has no thermal conductivity for element ", elem_id);

        const Precision length = length_reference();
        logging::error(length > Precision(0),
            "T3D2: zero reference length for element ", elem_id);

        const Precision scale = material->get_thermal_conductivity()
                              * get_section()->area_ / length;

        MapMatrix matrix(buffer, 2, 2);
        matrix << scale, -scale,
                 -scale,  scale;
        return matrix;
    }

    MapMatrix capacity(Precision* buffer) override {
        auto material = get_material();
        logging::error(material->has_density(),
            "T3D2: material has no density for element ", elem_id);
        logging::error(material->has_thermal_specific_heat(),
            "T3D2: material has no specific heat for element ", elem_id);

        const Precision length = length_reference();
        const Precision scale = material->get_density()
                              * material->get_thermal_specific_heat()
                              * get_section()->area_ * length / Precision(6);

        MapMatrix matrix(buffer, 2, 2);
        matrix << Precision(2) * scale, scale,
                  scale, Precision(2) * scale;
        return matrix;
    }

    void compute_heat_flux(Field& heat_flux, const Field& temperature) override {
        logging::error(heat_flux.domain == FieldDomain::ELEMENT_IP,
            "T3D2: heat flux output must use ELEMENT_IP domain");
        logging::error(heat_flux.components >= 3,
            "T3D2: heat flux output requires three components");
        logging::error(temperature.domain == FieldDomain::NODE && temperature.components == 1,
            "T3D2: temperature must be a scalar NODE field");

        auto material = get_material();
        logging::error(material->has_thermal_conductivity(),
            "T3D2: material has no thermal conductivity for element ", elem_id);

        const Precision length = length_reference();
        logging::error(length > Precision(0),
            "T3D2: zero reference length for element ", elem_id);

        const Precision t0 = temperature(static_cast<Index>(node_ids[0]), 0);
        const Precision t1 = temperature(static_cast<Index>(node_ids[1]), 0);
        const Vec3 flux = -material->get_thermal_conductivity()
                        * (t1 - t0) / length * direction_reference();

        const Index row = ip_index(0);
        for (Dim component = 0; component < 3; ++component) {
            heat_flux(row, component) = flux(component);
        }
    }
};

} // namespace fem::model
