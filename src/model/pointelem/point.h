//
// Created by f_eggers on 11.12.2024.
//

#ifndef POINT_H
#define POINT_H

#include "../element/element_structural.h"
#include "../../core/types_cls.h"

namespace fem::model {

struct Point : fem::model::StructuralElement {

    ID node_id;

    explicit Point(ID p_elem_id, ID p_node_id)
        : StructuralElement(p_elem_id)
        , node_id(p_node_id) {}

    std::string type_name() const override { return "P"; }

    ElDofs     dofs() override {
        return ElDofs {false, false, false, false, false, false};
    }
    Dim        dimensions() override {
        return 3;
    }
    Dim        n_nodes() override {
        return 1;
    };
    Dim        n_integration_points() override{
        return 0;
    };
    ID*        nodes() override{
        return &node_id;
    };
    SurfacePtr surface(ID surface_id) override {
        return nullptr;
    }
    Precision  volume() override {
        return 0;
    };
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        throw std::runtime_error("Not implemented");
    }

    ;
    MapMatrix mass(Precision* buffer) override;
    ;
    void compute_stress_strain_nodal(Field& displacement, Field& stress, Field& strain) override {
        (void) displacement;
        (void) stress;
        (void) strain;
    };
    void compute_stress_strain(Field& ip_stress, Field& ip_strain, Field& displacement, int ip_offset) override {
        (void) ip_stress;
        (void) ip_strain;
        (void) displacement;
        (void) ip_offset;
    };
    void apply_vload(Field& node_loads, Vec3 load) override {
        (void) node_loads;
        (void) load;
    };
    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    };
    void compute_compliance(Field& displacement, Field& result) override {
        (void) displacement;
        (void) result;
    };
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override {
        (void) displacement;
        (void) result;
    };
    Stresses stress(Field& displacement, std::vector<Vec3>& rst) override {
        return {};
    };
    Strains  strain(Field& displacement, std::vector<Vec3>& rst) override {
        return {};
    };
};

}



#endif //POINT_H
