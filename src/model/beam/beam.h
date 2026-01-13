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
#include <cmath>

namespace fem {
namespace model {

/**
 * @struct BeamElement
 * @brief Provides shared functionality for beam formulations with `N` nodes.
 */
template<Index N>
struct BeamElement : StructuralElement {
    std::array<ID, N> node_ids{}; ///< Connectivity of the beam element.
    ID orientation_node_id_ = static_cast<ID>(-1); ///< Optional node encoding the n1 direction.

    BeamElement(ID elem_id, std::array<ID, N> node_ids_in, ID orientation_node_id = static_cast<ID>(-1))
        : StructuralElement(elem_id)
        , node_ids(node_ids_in)
        , orientation_node_id_(orientation_node_id) {}

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

    ID orientation_node() const { return orientation_node_id_; }

    bool has_orientation_node() const { return orientation_node_id_ >= 0; }

    Vec3 orientation_direction() {
        constexpr Precision kOrientationEps = static_cast<Precision>(1e-12);
        BeamSection* section = get_section();
        const bool section_has_direction = section && section->n1.norm() > kOrientationEps;
        const bool element_has_node = has_orientation_node();

        logging::error(element_has_node || section_has_direction,
                       "Beam element ", this->elem_id,
                       " requires either an orientation node or a section-defined n1 vector");

        Vec3 n1_vec = Vec3::Zero();
        if (element_has_node) {
            auto row = this->_model_data->node_data.get(POSITION).row(orientation_node_id_);
            n1_vec = Vec3(row(0), row(1), row(2));
        } else {
            n1_vec = section->n1;
        }

        logging::error(n1_vec.norm() > kOrientationEps,
                       "Orientation vector for element ", this->elem_id, " must be non-zero");
        return n1_vec.normalized();
    }

    Precision length() {
        Precision l = 0;
        for (Index i = 0; i < N - 1; i++) {
            l += (coordinate(i) - coordinate(i + 1)).norm();
        }
        return l;
    }

    Precision volume() override { return get_profile()->A * length(); }

    // Base rotation (provided section frame; no principal rotation)
    Mat3 base_rotation_matrix() {
        Vec3 x = (coordinate(1) - coordinate(0)).normalized();
        Vec3 y = orientation_direction();
        Vec3 z = x.cross(y).normalized();
        y = z.cross(x).normalized();

        Mat3 mat{};
        mat(0, 0) = x(0); mat(0, 1) = x(1); mat(0, 2) = x(2);
        mat(1, 0) = y(0); mat(1, 1) = y(1); mat(1, 2) = y(2);
        mat(2, 0) = z(0); mat(2, 1) = z(1); mat(2, 2) = z(2);
        return mat;
    }

    // Compute principal rotation angle about local x from section inertias
    Precision principal_angle() {
        Profile* pr = get_profile();
        const Precision Iy  = pr->I_y;
        const Precision Iz  = pr->I_z;
        const Precision Iyz = pr->I_yz;

        const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
        if (std::abs(Iyz) <= scale * Precision(1e-14)) return Precision(0);
        return Precision(0.5) * std::atan2(Precision(2) * Iyz, Iy - Iz);
    }

    // Rotation aligning y/z to principal bending axes
    Mat3 principal_rotation_matrix() {
        Mat3 Rb = base_rotation_matrix();
        const Precision phi = principal_angle();
        if (phi == Precision(0)) return Rb;

        Eigen::Matrix<Precision, 1, 3> rx = Rb.row(0);
        Eigen::Matrix<Precision, 1, 3> ry = Rb.row(1);
        Eigen::Matrix<Precision, 1, 3> rz = Rb.row(2);

        const Precision c = std::cos(phi);
        const Precision s = std::sin(phi);

        Eigen::Matrix<Precision, 1, 3> ry_p =  c * ry + s * rz;
        Eigen::Matrix<Precision, 1, 3> rz_p = -s * ry + c * rz;

        Mat3 Rp = Rb;
        Rp.row(0) = rx;
        Rp.row(1) = ry_p;
        Rp.row(2) = rz_p;
        return Rp;
    }

    // Backward-compatible alias used by internal computations
    Mat3 rotation_matrix() { return principal_rotation_matrix(); }

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

    // Transformation using original section frame (no principal rotation)
    StaticMatrix<N * 6, N * 6> transformation_base() {
        StaticMatrix<N * 6, N * 6> T;
        T.setZero();
        Mat3 R = base_rotation_matrix();

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
        const Precision L = length();
        const Precision A = get_profile()->A;
        if (L <= Precision(0) || A <= Precision(0)) return;
        const Vec3 F = load * (A * L / Precision(N));
        for (Index i = 0; i < N; ++i) {
            const ID n_id = node_ids[i];
            node_loads(n_id, 0) += F(0);
            node_loads(n_id, 1) += F(1);
            node_loads(n_id, 2) += F(2);
        }
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

    std::vector<Vec6> section_forces(NodeData& displacement) override {
        std::vector<Vec6> result;
        result.resize(N);

        // 1) globale Verschiebungen einsammeln: u_global (N*6 x 1)
        Eigen::Matrix<Precision, N * 6, 1> u_global;
        for (Index i = 0; i < N; ++i) {
            const ID nid = node_ids[i];
            auto row = displacement.row(nid); // [ux, uy, uz, rx, ry, rz]
            for (Index d = 0; d < 6; ++d) {
                u_global(i * 6 + d) = row(d);
            }
        }

        // 2) globale Steifigkeit und Transformation holen
        const auto K_global = stiffness_impl(); // K_global = T^T * K_local * T (T uses principal frame)
        const auto T        = transformation_base(); // For output, use original (non-rotated) local frame

        // 3) globale Knotenkräfte
        const auto f_global = K_global * u_global; // N*6 x 1

        // 4) lokale Schnittgrößen (N, Vy, Vz, T, My, Mz) im Balkensystem
        const auto q_local = T * f_global; // N*6 x 1, in original local frame

        // 5) in N Vec6-Blöcke zuschneiden
        for (Index i = 0; i < N; ++i) {
            Vec6 q_i;
            for (Index d = 0; d < 6; ++d) {
                q_i(d) = q_local(i * 6 + d) * ((i == 1 && N==2) ? 1:-1);
            }
            result[i] = q_i;
        }

        return result;
    }


};

} // namespace model
} // namespace fem

