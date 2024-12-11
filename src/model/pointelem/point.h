//
// Created by f_eggers on 11.12.2024.
//

#ifndef POINT_H
#define POINT_H
#include "../element/element_structural.h"

namespace fem::model {

struct Point : fem::model::StructuralElement {

    ID node_id;

    explicit Point(ID p_elem_id, ID p_node_id)
        : StructuralElement(p_elem_id)
        , node_id(p_node_id) {}

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
    MapMatrix  stiffness(Precision* buffer) override {
        return MapMatrix(buffer, 6, 6).setZero();
    };
    MapMatrix  mass(Precision* buffer) override {
        logging::error(_section != nullptr, "Section not assigned to element ", elem_id);
        logging::error(_section->as<PointMassSection>(), "Section is not a PointMassSection");
        PointMassSection* section = _section->as<PointMassSection>();
        MapMatrix res = MapMatrix(buffer, 6, 6).setZero();
        res(0, 0) = section->mass;
        res(1, 1) = section->mass;
        res(2, 2) = section->mass;
        res(3, 3) = section->rotary_inertia(0);
        res(4, 4) = section->rotary_inertia(1);
        res(5, 5) = section->rotary_inertia(2);
        return res;
    };
    void compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) override {
        (void) displacement;
        (void) stress;
        (void) strain;
    };
    void compute_stress_strain(NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) override {
        (void) displacement;
        (void) stress;
        (void) strain;
        (void) xyz;
    };
    void apply_vload(NodeData& node_loads, Vec3 load) override {
        (void) node_loads;
        (void) load;
    };
    void apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    };
    void compute_compliance(NodeData& displacement, ElementData& result) override {
        (void) displacement;
        (void) result;
    };
    void compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) override {
        (void) displacement;
        (void) result;
    };
};

}



#endif //POINT_H
