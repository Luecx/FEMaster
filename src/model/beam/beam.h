/**
 * @file beam.h
 * @brief Declares the templated base class for structural beam elements.
 *
 * `BeamElement<N>` factors out common geometry and section/material handling for
 * concrete beam implementations while deriving from `StructuralElement`.
 *
 * @see src/model/beam/b33.h
 */

#pragma once

#include "../../core/core.h"
#include "../../section/section_beam.h"
#include "../element/element_structural.h"
#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"

#include <array>

namespace fem {
namespace model {

/**
 * @struct BeamElement
 * @brief Provides shared functionality for beam formulations with `N` nodes.
 */
template<Index N>
struct BeamElement : StructuralElement {
    std::array<ID, N> node_ids{}; ///< Connectivity of the beam element.

    BeamElement(ID elem_id, std::array<ID, N> node_ids_in)
        : StructuralElement(elem_id)
        , node_ids(node_ids_in) {}

    ~BeamElement() override = default;

    BeamSection* get_section() {
        if (!this->_section) {
            logging::error(false, "Section not set for element ", this->elem_id);
        }
        if (!this->_section->template as<BeamSection>()) {
            logging::error(false, "Section is not a beam section for element ", this->elem_id);
        }
        return this->_section->template as<BeamSection>();
    }

    Profile* get_profile() { return get_section()->profile.get(); }

    material::MaterialPtr get_material() {
        BeamSection* section = get_section();
        if (!section->material) {
            logging::error(false, "Material not set for element ", this->elem_id);
        }
        return section->material;
    }

    material::IsotropicElasticity* get_elasticity() {
        BeamSection* section = get_section();
        if (!section->material) {
            logging::error(false, "Material not set for element ", this->elem_id);
        }
        if (!section->material->has_elasticity()) {
            logging::error(false, "Material has no elasticity assigned");
        }
        if (!section->material->elasticity()->template as<material::IsotropicElasticity>()) {
            logging::error(false, "Material is not isotropic for element ", this->elem_id);
        }
        return section->material->elasticity()->template as<material::IsotropicElasticity>();
    }

    Vec3 coordinate(Index index) {
        auto node_id = node_ids[index];
        auto row = this->_model_data->node_data.get(POSITION).row(node_id);
        return Vec3(row(0), row(1), row(2));
    }

    Precision length() {
        Precision l = 0;
        for (Index i = 0; i < N - 1; i++) {
            l += (coordinate(i) - coordinate(i + 1)).norm();
        }
        return l;
    }

    Precision volume() override { return get_profile()->A * length(); }

    Mat3 rotation_matrix() {
        Vec3 x = (coordinate(1) - coordinate(0)).normalized();
        Vec3 y = get_section()->n1;
        Vec3 z = x.cross(y).normalized();
        y = z.cross(x).normalized();

        Mat3 mat{};
        mat(0, 0) = x(0);
        mat(0, 1) = x(1);
        mat(0, 2) = x(2);
        mat(1, 0) = y(0);
        mat(1, 1) = y(1);
        mat(1, 2) = y(2);
        mat(2, 0) = z(0);
        mat(2, 1) = z(1);
        mat(2, 2) = z(2);
        return mat;
    }

    virtual StaticMatrix<N * 6, N * 6> stiffness_impl() = 0;
    virtual StaticMatrix<N * 6, N * 6> stiffness_geom_impl(IPData& ip_stress, int offset) = 0;
    virtual StaticMatrix<N * 6, N * 6> mass_impl() = 0;

    StaticMatrix<N * 6, N * 6> transformation() {
        StaticMatrix<N * 6, N * 6> T;
        T.setZero();
        Mat3 R = rotation_matrix();

        for (Index i = 0; i < N; i++) {
            for (Dim j = 0; j < 3; j++) {
                for (Dim k = 0; k < 3; k++) {
                    T(i * 6 + j, i * 6 + k) = R(j, k);
                    T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
                }
            }
        }
        return T;
    }

    MapMatrix stiffness(Precision* buffer) override {
        MapMatrix result(buffer, N * 6, N * 6);
        result = stiffness_impl();
        return result;
    }

    MapMatrix stiffness_geom(Precision* buffer, IPData& ip_stress, int ip_start_idx) override {
        MapMatrix result(buffer, N * 6, N * 6);
        result = stiffness_geom_impl(ip_stress, ip_start_idx);
        return result;
    }

    MapMatrix mass(Precision* buffer) override {
        MapMatrix result(buffer, N * 6, N * 6);
        result = mass_impl();
        return result;
    }

    void compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) override {
        (void)displacement;
        (void)stress;
        (void)strain;
    }

    void compute_stress_strain(IPData& ip_stress, IPData& ip_strain, NodeData& displacement, int ip_offset) override {
        (void)ip_stress;
        (void)ip_strain;
        (void)displacement;
        (void)ip_offset;
    }

    void apply_vload(NodeData& node_loads, Vec3 load) override {
        (void)node_loads;
        (void)load;
    }

    void apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) override {
        (void)node_loads;
        (void)node_temp;
        (void)ref_temp;
    }

    void compute_compliance(NodeData& displacement, ElementData& result) override {
        (void)displacement;
        (void)result;
    }

    void compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) override {
        (void)displacement;
        (void)result;
    }

    ElDofs dofs() override { return ElDofs{true, true, true, true, true, true}; }
    Dim dimensions() override { return 3; }
    Dim n_nodes() override { return N; }
    Dim n_integration_points() override { return 1; }
    ID* nodes() override { return node_ids.data(); }
    SurfacePtr surface(ID surface_id) override {
        (void)surface_id;
        return nullptr;
    }
    Stresses stress(NodeData& displacement, std::vector<Vec3>& rst) override {
        (void)displacement;
        (void)rst;
        return {};
    }
    Strains strain(NodeData& displacement, std::vector<Vec3>& rst) override {
        (void)displacement;
        (void)rst;
        return {};
    }
};

} // namespace model
} // namespace fem

