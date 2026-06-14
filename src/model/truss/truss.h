#pragma once

#include "../../core/core.h"
#include "../element/element_structural.h"
#include "../../section/section_truss.h"
#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"

namespace fem::model {
template<Index N>
struct TrussElement : StructuralElement {
    static_assert(N == 2, "Only 2-node truss elements are supported.");

    std::array<ID, N> node_ids;

    TrussElement(ID p_elem_id, std::array<ID, N> p_node_ids)
        : StructuralElement(p_elem_id), node_ids{p_node_ids} {}

    std::string type_name() const override { return (N == 2 ? std::string("T3") : std::string("TRUSS")); }

    ~TrussElement() override = default;

    // --- Section and material access ---
    TrussSection* get_section() {
        if (!_section || !_section->template as<TrussSection>())
            logging::error(false, "TrussElement: invalid or missing TrussSection for element ", elem_id);
        return _section->template as<TrussSection>();
    }

    material::MaterialPtr get_material() {
        auto section = get_section();
        if (!section->material_)
            logging::error(false, "TrussElement: no material set for element ", elem_id);
        return section->material_;
    }

    material::IsotropicElasticity* get_elasticity() {
        auto mat = get_material();
        if (!mat->has_elasticity())
            logging::error(false, "TrussElement: material has no elasticity.");
        return mat->elasticity()->template as<material::IsotropicElasticity>();
    }

    // --- Geometry ---
    Vec3 coordinate(Index index) {
        return this->node_position(node_ids[index]);
    }

    Precision length() {
        return (coordinate(1) - coordinate(0)).norm();
    }

    Vec3 direction() {
        return (coordinate(1) - coordinate(0)).normalized();
    }

    Precision volume() override {
        return get_section()->area_ * length();
    }

    // --- Stiffness ---
    virtual StaticMatrix<N * 3, N * 3> stiffness_impl() {
        Precision E = get_elasticity()->youngs;
        Precision A = get_section()->area_;
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

        const auto mat = get_material();
        if (!mat->has_density()) return M;

        Precision rho = mat->get_density();
        Precision A = get_section()->area_;
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
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        logging::error(ip_stress.components >= 1,
                       "TrussElement: geometric stiffness requires nonlinear IP stress component 0");

        const Precision L = length();
        logging::error(L > Precision(0),
                       "TrussElement: zero current length in stiffness_geom for element ",
                       this->elem_id);

        const Vec3 n = direction();
        const Precision sigma = ip_stress(static_cast<Index>(ip_start_idx), 0);
        const Precision axial_force = get_section()->area_ * sigma;
        const Mat3 projector = Mat3::Identity() - n * n.transpose();
        const Mat3 block = (axial_force / L) * projector;

        StaticMatrix<N * 3, N * 3> Kg = StaticMatrix<N * 3, N * 3>::Zero();
        Kg.template block<3, 3>(0, 0) = block;
        Kg.template block<3, 3>(0, 3) = -block;
        Kg.template block<3, 3>(3, 0) = -block;
        Kg.template block<3, 3>(3, 3) = block;

        MapMatrix result(buffer, N * 3, N * 3);
        result = Kg;
        return result;
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

    RowMatrix stress_strain_nodal_rst() override {
        RowMatrix rst(N, 3);
        rst.setZero();
        rst(0, 0) = Precision(-1);
        rst(1, 0) = Precision(1);
        return rst;
    }

    RowMatrix stress_strain_ip_rst() override {
        RowMatrix rst(1, 3);
        rst.setZero();
        return rst;
    }

    void compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        logging::error(strain != nullptr || stress != nullptr,
                       "TrussElement: compute_stress_strain requires at least one output field");
        logging::error(rst.cols() >= 1,
                       "TrussElement: stress/strain evaluation coordinates require at least 1 column");

        Precision strain_value = Precision(0);
        Precision stress_value = Precision(0);

        if (use_green_lagrange_nl) {
            const Vec3 X0 = coordinate(0);
            const Vec3 X1 = coordinate(1);
            const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
            const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
            const Vec3 x0 = X0 + u0;
            const Vec3 x1 = X1 + u1;

            const Precision L0 = (X1 - X0).norm();
            const Precision l = (x1 - x0).norm();
            logging::error(L0 > Precision(0),
                           "TrussElement: zero reference length in compute_stress_strain for element ",
                           this->elem_id);
            logging::error(l > Precision(0),
                           "TrussElement: zero current length in compute_stress_strain for element ",
                           this->elem_id);

            const Precision lambda = l / L0;
            strain_value = Precision(0.5) * (lambda * lambda - Precision(1));
            const Precision second_piola_stress = get_elasticity()->youngs * strain_value;
            stress_value = lambda * second_piola_stress;
        } else {
            const Precision L = length();
            logging::error(L > Precision(0), "TrussElement: zero length in compute_stress_strain for element ",
                           this->elem_id);
            const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
            const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
            strain_value = (u1 - u0).dot(direction()) / L;
            stress_value = get_elasticity()->youngs * strain_value;
        }

        for (Index i = 0; i < rst.rows(); ++i) {
            const Index row = static_cast<Index>(offset) + i;
            if (strain) {
                for (Index j = 0; j < strain->components; ++j) (*strain)(row, j) = Precision(0);
                (*strain)(row, 0) = strain_value;
            }
            if (stress) {
                for (Index j = 0; j < stress->components; ++j) (*stress)(row, j) = Precision(0);
                (*stress)(row, 0) = stress_value;
            }
        }
    }

    void compute_internal_force_nonlinear(Field& node_forces,
                                          const Field& ip_stress,
                                          int ip_offset) override {
        logging::error(node_forces.domain == FieldDomain::NODE,
                       "TrussElement: nonlinear internal force output must use NODE domain");
        logging::error(node_forces.components >= 3,
                       "TrussElement: nonlinear internal force output requires at least 3 components");
        logging::error(ip_stress.components >= 1,
                       "TrussElement: nonlinear internal force requires IP stress component 0");

        const Precision L = length();
        logging::error(L > Precision(0),
                       "TrussElement: zero current length in compute_internal_force_nonlinear for element ",
                       this->elem_id);

        const Vec3 n = direction();
        const Precision sigma = ip_stress(static_cast<Index>(ip_offset), 0);
        const Vec3 force = get_section()->area_ * sigma * n;

        const Index n0 = static_cast<Index>(node_ids[0]);
        const Index n1 = static_cast<Index>(node_ids[1]);
        for (Dim d = 0; d < 3; ++d) {
            node_forces(n0, d) -= force(d);
            node_forces(n1, d) += force(d);
        }
    }

    bool compute_beam_section_forces(Field& section_forces,
                                     const Field& displacement,
                                     int offset) override {
        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

        const Precision L = length();
        logging::error(L > Precision(0),
                       "TrussElement: zero length in compute_beam_section_forces for element ",
                       this->elem_id);

        const Vec3 t = direction();
        const Precision axial_strain = (u1 - u0).dot(t) / L;
        const Precision axial_force = get_elasticity()->youngs * get_section()->area_ * axial_strain;

        for (Index i = 0; i < N; ++i) {
            for (Index d = 0; d < section_forces.components; ++d) {
                section_forces(static_cast<Index>(offset) + i, d) = Precision(0);
            }
            section_forces(static_cast<Index>(offset) + i, 0) = axial_force;
        }
        return true;
    }
    void apply_vload(Field&, Vec3) override {}
    void apply_tload(Field&, const Field&, Precision) override {}
    void compute_compliance(Field&, Field&) override {}
    void compute_compliance_angle_derivative(Field&, Field&) override {}

    // Integrate vector field via 1-point midpoint rule and equal nodal distribution
    void integrate_vec_field(Field& node_loads,
                             bool scale_by_density,
                             const VecField& field) override {
        const Precision L = length();
        const Precision A = get_section()->area_;
        if (L <= Precision(0) || A <= Precision(0)) return;

        // Midpoint
        Vec3 x_mid = (coordinate(0) + coordinate(1)) * Precision(0.5);

        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "TrussElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }

        const Vec3 F = field(x_mid) * (rho * A * L);
        for (Index i = 0; i < N; ++i) {
            const ID n_id = node_ids[i];
            node_loads(n_id, 0) += F(0) * Precision(0.5);
            node_loads(n_id, 1) += F(1) * Precision(0.5);
            node_loads(n_id, 2) += F(2) * Precision(0.5);
        }
    }

    // --- Generic integrations (midpoint over volume A*L) ---
    Precision integrate_scalar_field(bool scale_by_density,
                                     const ScalarField& field) override {
        const Precision L = length();
        const Precision A = get_section()->area_;
        if (L <= Precision(0) || A <= Precision(0)) return Precision(0);
        Vec3 x_mid = (coordinate(0) + coordinate(1)) * Precision(0.5);
        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "TrussElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }

    Vec3 integrate_vector_field(bool scale_by_density,
                                const VecField& field) override {
        const Precision L = length();
        const Precision A = get_section()->area_;
        if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();
        Vec3 x_mid = (coordinate(0) + coordinate(1)) * Precision(0.5);
        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "TrussElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }

    Mat3 integrate_tensor_field(bool scale_by_density,
                                const TenField& field) override {
        const Precision L = length();
        const Precision A = get_section()->area_;
        if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();
        Vec3 x_mid = (coordinate(0) + coordinate(1)) * Precision(0.5);
        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "TrussElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }
};

using T3 = TrussElement<2>;
}  // namespace fem::model
