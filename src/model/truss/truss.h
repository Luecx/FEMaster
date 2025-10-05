#pragma once

#include "../../core/core.h"
#include "../element/element_structural.h"
#include "../../section/section_beam.h"
#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"

namespace fem::model {

template<Index N>
struct TrussElement : StructuralElement {

    static_assert(N == 2, "Only 2-node truss elements are supported.");

    std::array<ID, N> node_ids;

    TrussElement(ID p_elem_id, std::array<ID, N> p_node_ids)
        : StructuralElement(p_elem_id), node_ids{p_node_ids} {}

    ~TrussElement() override = default;

    // --- Section and material access ---
    BeamSection* get_section() {
        if (!_section || !_section->template as<BeamSection>())
            logging::error(false, "TrussElement: invalid or missing BeamSection for element ", elem_id);
        return _section->template as<BeamSection>();
    }

    Profile* get_profile() {
        return get_section()->profile.get();
    }

    material::MaterialPtr get_material() {
        auto section = get_section();
        if (!section->material)
            logging::error(false, "TrussElement: no material set for element ", elem_id);
        return section->material;
    }

    material::IsotropicElasticity* get_elasticity() {
        auto mat = get_material();
        if (!mat->has_elasticity())
            logging::error(false, "TrussElement: material has no elasticity.");
        return mat->elasticity()->template as<material::IsotropicElasticity>();
    }

    // --- Geometry ---
    Vec3 coordinate(Index index) {
        auto node_id = node_ids[index];
        auto row = this->_model_data->node_data.get(POSITION).row(node_id);
        return Vec3(row(0), row(1), row(2));
    }

    Precision length() {
        return (coordinate(1) - coordinate(0)).norm();
    }

    Vec3 direction() {
        return (coordinate(1) - coordinate(0)).normalized();
    }

    Precision volume() override {
        return get_profile()->A * length();
    }

    // --- Stiffness ---
    virtual StaticMatrix<N * 3, N * 3> stiffness_impl() {
        Precision E = get_elasticity()->youngs;
        Precision A = get_profile()->A;
        Precision L = length();

        StaticMatrix<2, 2> k_local;
        k_local << 1, -1,
                  -1,  1;

        Vec3 t = direction();

        StaticMatrix<3 * N, 2> P = StaticMatrix<3 * N, 2>::Zero();
        P.block(0, 0, 3, 1) = t;
        P.block(3, 1, 3, 1) = t;

        return P * (k_local * (E * A / L)) * P.transpose();
    }

    // --- Mass ---
    virtual StaticMatrix<N * 3, N * 3> mass_impl() {
        StaticMatrix<N * 3, N * 3> M = StaticMatrix<N * 3, N * 3>::Zero();

        Precision rho = get_material()->get_density();
        Precision A = get_profile()->A;
        Precision L = length();
        Precision m = rho * A * L;

        for (Index i = 0; i < N; ++i)
            M.block(i * 3, i * 3, 3, 3) = Mat3::Identity() * (m / 2.0);

        return M;
    }

    // --- Interface overrides ---
    MapMatrix stiffness(Precision* buffer) override {
        MapMatrix result(buffer, N * 3, N * 3);
        result = stiffness_impl();
        return result;
    }

    // --- Interface overrides ---
    MapMatrix stiffness_geom(Precision* buffer, IPData& ip_stress, int ip_start_idx) override {
        throw std::runtime_error("Not implemented yet");
    }

    MapMatrix mass(Precision* buffer) override {
        MapMatrix result(buffer, N * 3, N * 3);
        result = mass_impl();
        return result;
    }

    ElDofs dofs() override {
        return ElDofs{true, true, true, false, false, false};
    }

    Dim dimensions() override {
        return 3;
    }

    Dim n_nodes() override {
        return N;
    }

    Dim n_integration_points() override {
        return 1;
    }

    ID* nodes() override {
        return node_ids.data();
    }

    SurfacePtr surface(ID) override {
        return nullptr;
    }

    // --- Stress/Strain ---
    Stresses stress(NodeData& displacement, std::vector<Vec3>& rst) override {
        throw std::runtime_error("Not implemented yet");
        // Vec3 u0 = displacement.row(node_ids[0]);
        // Vec3 u1 = displacement.row(node_ids[1]);
        //
        // Vec3 t = direction();
        // Precision L = length();
        //
        // Precision axial_strain = (u1 - u0).dot(t) / L;
        // Precision stress_val = get_elasticity()->youngs * axial_strain;
        //
        // if (!rst.empty()) rst[0] = 0.5 * (coordinate(0) + coordinate(1));
        //
        // Stresses result(1, 1);
        // result(0, 0) = stress_val;
        // return result;
    }

    Strains strain(NodeData& displacement, std::vector<Vec3>& rst) override {
        throw std::runtime_error("Not implemented yet");

        // Vec3 u0 = displacement.row(node_ids[0]);
        // Vec3 u1 = displacement.row(node_ids[1]);
        //
        // Vec3 t = direction();
        // Precision L = length();
        //
        // Precision axial_strain = (u1 - u0).dot(t) / L;
        //
        // if (!rst.empty()) rst[0] = 0.5 * (coordinate(0) + coordinate(1));
        //
        // Strains result(1, 1);
        // result(0, 0) = axial_strain;
        // return result;
    }

    void compute_stress_strain(IPData& ip_stress, IPData& ip_strain, NodeData& displacement, int ip_offset) override {
        (void) ip_stress;
        (void) ip_strain;
        (void) displacement;
        (void) ip_offset;
    }
    void compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) override {
        (void) displacement;
        (void) stress;
        (void) strain;
    }
    void apply_vload(NodeData&, Vec3) override {}
    void apply_tload(NodeData&, NodeData&, Precision) override {}
    void compute_compliance(NodeData&, ElementData&) override {}
    void compute_compliance_angle_derivative(NodeData&, ElementData&) override {}
};

using T3 = TrussElement<2>;

}  // namespace fem::model
