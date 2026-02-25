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

    Vec3 coordinate(Index index) { return this->node_position(node_ids[index]); }

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
            n1_vec = this->node_position(orientation_node_id_);
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

    // --- Kinematic rigid-offset (local y/z offsets) -------------------------

    /**
     * @brief Builds the 6x6 kinematic rigid-offset matrix B for a point shift d=(0,dy,dz) in local coordinates.
     *
     * For u = [ux,uy,uz, rx,ry,rz]^T:
     *   u_A = B(dy,dz) * u_R  with  u_A = u_R + theta x d
     */
    static StaticMatrix<6, 6> rigid_offset_6(Precision dy, Precision dz) {
        StaticMatrix<6, 6> B = StaticMatrix<6, 6>::Identity();

        // u += theta x d, with d=(0,dy,dz)
        // [ux] +=  ry*dz - rz*dy
        // [uy] += -rx*dz
        // [uz] +=  rx*dy
        B(0, 4) =  dz;
        B(0, 5) = -dy;

        B(1, 3) = -dz;
        B(2, 3) =  dy;

        return B;
    }

    /**
     * @brief Builds the Nx6 by Nx6 block-diagonal rigid-offset matrix for all element nodes.
     */
    static StaticMatrix<N * 6, N * 6> rigid_offset_N(Precision dy, Precision dz) {
        StaticMatrix<N * 6, N * 6> B = StaticMatrix<N * 6, N * 6>::Identity();
        const StaticMatrix<6, 6> Bi = rigid_offset_6(dy, dz);
        for (Index a = 0; a < N; ++a) {
            B.template block<6, 6>(a * 6, a * 6) = Bi;
        }
        return B;
    }

    /**
     * @brief Rotate a (y,z) offset expressed in the base section frame into the principal frame used by rotation_matrix().
     *
     * This uses the same angle phi as the principal axis rotation:
     *   y_p =  c*y + s*z
     *   z_p = -s*y + c*z
     */
    static void rotate_yz_to_principal(Precision phi, Precision& y, Precision& z) {
        if (phi == Precision(0)) return;
        const Precision c = std::cos(phi);
        const Precision s = std::sin(phi);
        const Precision y_p =  c * y + s * z;
        const Precision z_p = -s * y + c * z;
        y = y_p;
        z = z_p;
    }

    // --- FEM interface ------------------------------------------------------

    virtual StaticMatrix<N * 6, N * 6> stiffness_impl() = 0;
    virtual StaticMatrix<N * 6, N * 6> stiffness_geom_impl(const Field& ip_stress, int offset) = 0;
    virtual StaticMatrix<N * 6, N * 6> mass_impl() = 0;

    StaticMatrix<N * 6, N * 6> transformation() {
        StaticMatrix<N * 6, N * 6> T;
        T.setZero();
        Mat3 R = rotation_matrix();

        for (Index i = 0; i < N; i++) {
            for (Dim j = 0; j < 3; j++) {
                for (Dim k = 0; k < 3; k++) {
                    T(i * 6 + j,     i * 6 + k)     = R(j, k);
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
                    T(i * 6 + j,     i * 6 + k)     = R(j, k);
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

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        MapMatrix result(buffer, N * 6, N * 6);
        result = stiffness_geom_impl(ip_stress, ip_start_idx);
        return result;
    }

    MapMatrix mass(Precision* buffer) override {
        MapMatrix result(buffer, N * 6, N * 6);
        result = mass_impl();
        return result;
    }

    void compute_stress_strain_nodal(Field& displacement, Field& stress, Field& strain) override {
        (void)displacement;
        (void)stress;
        (void)strain;
    }

    void compute_stress_strain(Field& ip_stress, Field& ip_strain, Field& displacement, int ip_offset) override {
        (void)ip_stress;
        (void)ip_strain;
        (void)displacement;
        (void)ip_offset;
    }

    void apply_vload(Field& node_loads, Vec3 load) override {
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

    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override {
        (void)node_loads;
        (void)node_temp;
        (void)ref_temp;
    }

    void compute_compliance(Field& displacement, Field& result) override {
        (void)displacement;
        (void)result;
    }

    void compute_compliance_angle_derivative(Field& displacement, Field& result) override {
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
    Stresses stress(Field& displacement, std::vector<Vec3>& rst) override {
        (void)displacement;
        (void)rst;
        return {};
    }
    Strains strain(Field& displacement, std::vector<Vec3>& rst) override {
        (void)displacement;
        (void)rst;
        return {};
    }

    std::vector<Vec6> section_forces(Field& displacement) override {
        std::vector<Vec6> result;
        result.resize(N);

        Eigen::Matrix<Precision, N * 6, 1> u_global;
        for (Index i = 0; i < N; ++i) {
            const ID nid = node_ids[i];
            const Vec6 row = displacement.row_vec6(static_cast<Index>(nid)); // [ux, uy, uz, rx, ry, rz]
            for (Index d = 0; d < 6; ++d) {
                u_global(i * 6 + d) = row(d);
            }
        }

        // global stiffness and local frame for output
        const auto K_global = stiffness_impl();
        const auto T_out    = transformation_base(); // output in base section frame

        const auto f_global = K_global * u_global;
        const auto q_local  = T_out * f_global; // nodal resultants in base local frame, about the element DOF reference line

        for (Index i = 0; i < N; ++i) {
            Vec6 q_i;
            for (Index d = 0; d < 6; ++d) {
                q_i(d) = q_local(i * 6 + d) * ((i == 1 && N == 2) ? 1 : -1);
            }
            result[i] = q_i;
        }

        return result;
    }

    // Integrate vector field via simple 1-point mid-span sampling and equal distribution
    void integrate_vec_field(Field& node_loads,
                             bool scale_by_density,
                             const VecField& field) override {
        const Precision L = length();
        const Precision A = get_profile()->A;
        if (L <= Precision(0) || A <= Precision(0)) return;

        Vec3 x_mid = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_mid += coordinate(i);
        x_mid /= static_cast<Precision>(N);

        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }

        const Vec3 F = field(x_mid) * (rho * A * L);
        const Precision share = Precision(1) / static_cast<Precision>(N);
        for (Index i = 0; i < N; ++i) {
            const ID n_id = node_ids[i];
            node_loads(n_id, 0) += share * F(0);
            node_loads(n_id, 1) += share * F(1);
            node_loads(n_id, 2) += share * F(2);
        }
    }

    Precision integrate_scalar_field(bool scale_by_density,
                                     const ScalarField& field) override {
        const Precision L = length();
        const Precision A = get_profile()->A;
        if (L <= Precision(0) || A <= Precision(0)) return Precision(0);

        Vec3 x_mid = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_mid += coordinate(i);
        x_mid /= static_cast<Precision>(N);

        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }

    Vec3 integrate_vector_field(bool scale_by_density,
                                const VecField& field) override {
        const Precision L = length();
        const Precision A = get_profile()->A;
        if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();

        Vec3 x_mid = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_mid += coordinate(i);
        x_mid /= static_cast<Precision>(N);

        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }

    Mat3 integrate_tensor_field(bool scale_by_density,
                                const TenField& field) override {
        const Precision L = length();
        const Precision A = get_profile()->A;
        if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();

        Vec3 x_mid = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_mid += coordinate(i);
        x_mid /= static_cast<Precision>(N);

        Precision rho = 1.0;
        if (scale_by_density) {
            auto mat = get_material();
            logging::error(mat && mat->has_density(),
                           "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
            rho = mat->get_density();
        }
        return field(x_mid) * (rho * A * L);
    }
};

} // namespace model
} // namespace fem