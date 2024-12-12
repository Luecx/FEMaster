//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_H
#define SHELL_H

#include "../../core/core.h"
#include "../element/element_structural.h"
#include "../../section/section_shell.h"

#include <memory>

namespace fem::model {
template<Index N>
struct ShellElement : StructuralElement {
    std::array<ID, N> node_ids;

    ShellElement(ID p_elem_id, std::array<ID, N> p_node_ids) : StructuralElement(p_elem_id), node_ids {p_node_ids} {}
    virtual ~ShellElement() override = default;

    ShellSection* get_section() {
        if (!this->_section) {
            logging::error(false, "Section not set for element ", this->elem_id);
        }
        if (!this->_section->as<ShellSection>()) {
            logging::error(false, "Section is not a beam section for element ", this->elem_id);
        }
        return this->_section->as<ShellSection>();
    }
    material::MaterialPtr get_material() {
        BeamSection* section = get_section();
        if (!section->material) {
            logging::error(false, "Material not set for element ", this->elem_id);
        }
        return section->material;
    }
    material::IsotropicElasticity* get_elasticity() {
        auto mat_ptr = get_material();
        if (!mat_ptr->has_elasticity()) {
            logging::error(false, "Material has no elasticity assigned");
        }
        if (!mat_ptr->elasticity()->template as<material::IsotropicElasticity>()) {
            logging::error(false, "Material is not isotropic for element ", this->elem_id);
        }
        return mat_ptr->elasticity()->template as<material::IsotropicElasticity>();
    }

    // left out for childs
    virtual SurfacePtr surface(ID surface_id) override = 0;
    virtual Precision  volume() override = 0;
    virtual MapMatrix  stiffness(Precision* buffer) override = 0;
    virtual MapMatrix  mass(Precision* buffer) override = 0;

    virtual const fem::quadrature::Quadrature& integration_scheme() const = 0;

    virtual Mat3 covariant(const Vec3& pos);
    virtual Mat3 contravariant(const Vec3& pos);



    //
    StaticMatrix<N, 3> node_coords_global() {
        auto node_coords_system = this->_model_data->get(POSITION);
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            for (Index j = 0; j < 3; j++) {
                res(i, j) = node_coords_system(this->node_ids[i], j);
            }
        }
        return res;
    }

    Dim        n_integration_points() {
        return integration_scheme().count();
    }

    ElDofs     dofs() override {
        return ElDofs{true, true, true, true, true, true};
    }
    Dim        dimensions() override {
        return 3;
    }
    Dim        n_nodes() override {
        return N;
    }
    ID*        nodes() override {
        return node_ids.data();
    }

    // not implemented
    void compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) override {};
    void compute_stress_strain(NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) override{};
    void apply_vload(NodeData& node_loads, Vec3 load) override{};
    void apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) override{};
    void compute_compliance(NodeData& displacement, ElementData& result) override{};
    void compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) override{};
};
}


#endif //SHELL_H
