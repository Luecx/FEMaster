/**
 * @file c3d8i.cpp
 * @brief Implements the incompatible-mode eight-node hexahedral solid.
 *
 * The formulation uses the full C3D8 quadrature together with thirteen local
 * enhanced deformation-gradient modes. The first nine modes are the three
 * vector-valued linear natural-coordinate modes; the remaining four are
 * volumetric rs, rt, st and rst modes. The Jacobian-ratio transformation makes
 * every enhanced mode satisfy the element patch-test orthogonality condition.
 *
 * Linear mechanics assembles the coupled nodal/enhanced block system and
 * eliminates the local parameters by static condensation. Geometrically
 * nonlinear mechanics applies the enhanced modes multiplicatively in the
 * reference configuration, solves the thirteen local stationarity equations by
 * Newton iteration and condenses the complete material and geometric tangent.
 *
 * The finite-strain kinematics use
 *
 *     F_bar = F_compatible (I + sum_m alpha_m H_m),
 *
 * so a superposed spatial rigid rotation acts as
 * F_bar -> Q F_bar without changing the local enhanced mode space.
 *
 * Thermal loading, stress recovery and geometric prestress stiffness reconstruct
 * the same stationary enhanced state. Material history remains state-neutral in
 * all auxiliary paths; only the physical nonlinear tangent writes trial state.
 *
 * @see C3D8I
 * @see C3D8
 * @see SolidElement
 *
 * @author Finn Eggers
 * @date 02.10.2026
 */

#include "c3d8i.h"

#include "../../cos/rectangular_system.h"

#include <Eigen/LU>

#include <cmath>

namespace fem::model {

C3D8I::C3D8I(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8I::type_name() const {
    return "C3D8I";
}

/**
 * Builds the thirteen incompatible deformation-gradient basis tensors.
 *
 * Let J(xi) be the reference isoparametric Jacobian and J0 its value at the
 * element center. The nine principal modes are transformed from the natural
 * coordinate gradients
 *
 *     alpha_i * xi_i
 *
 * with the factor det(J0)/det(J). Four additional scalar modes
 * rs, rt, st and rst multiply the reference identity tensor. The determinant
 * ratio makes the physical-volume mean of every enhanced mode vanish.
 *
 * In the finite-strain path these tensors form the right multiplicative
 * enhancement of the compatible deformation gradient. They therefore remain in
 * the reference configuration and do not rotate independently with the current
 * spatial frame.
 *
 * @param reference_coords Global nodal coordinates in the reference configuration.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Thirteen reference-configuration enhanced gradient tensors.
 */
C3D8I::EnhancedModes C3D8I::enhanced_gradient_modes(
    const StaticMatrix<N, D>& reference_coords,
    Precision                 r,
    Precision                 s,
    Precision                 t
) {
    const Mat3      J  = this->jacobian(reference_coords, r, s, t);
    const Mat3      J0 = this->jacobian(reference_coords, Precision(0), Precision(0), Precision(0));
    const Precision j  = J.determinant();
    const Precision j0 = J0.determinant();

    logging::error(std::isfinite(j) && j > Precision(0),
        "C3D8I: invalid reference determinant in element ", elem_id,
        "\ndet(J): ", j);
    logging::error(std::isfinite(j0) && j0 > Precision(0),
        "C3D8I: invalid center reference determinant in element ", elem_id,
        "\ndet(J0): ", j0);

    const Precision                scale                = j0 / j;
    const Mat3                     natural_to_reference = J0.inverse().transpose();
    const std::array<Precision, 3> xi                   {r, s, t};

    EnhancedModes modes;
    for (auto& mode : modes) {
        mode.setZero();
    }

    // Nine principal modes: three displacement-vector components for each
    // natural-coordinate direction. Each tensor has only one nonzero row,
    // equal to the scaled reference dual vector of that natural direction.
    for (Index direction = 0; direction < 3; ++direction) {
        for (Index component = 0; component < 3; ++component) {
            modes[3 * direction + component].row(component) =
                (scale * xi[direction]) * natural_to_reference.row(direction);
        }
    }

    // Four volumetric modes used by C3D8I to improve the nearly incompressible
    // response without introducing any global pressure degree of freedom.
    const std::array<Precision, 4> theta {
        r * s,
        r * t,
        s * t,
        r * s * t
    };

    for (Index mode = 0; mode < 4; ++mode) {
        modes[9 + mode] = scale * theta[mode] * Mat3::Identity();
    }

    return modes;
}

/**
 * Converts the enhanced gradient tensors into infinitesimal engineering-strain
 * columns.
 *
 * Each column is the symmetric part of one enhanced displacement gradient in
 * the FEMaster volume-strain ordering
 * [eps_xx, eps_yy, eps_zz, gamma_yz, gamma_xz, gamma_xy]. The matrix is used by
 * linear stiffness, thermal loading, prestress recovery and linear output.
 *
 * @param modes Enhanced gradient tensors in the reference configuration.
 * @return Engineering-strain matrix with one column per enhanced parameter.
 */
C3D8I::Matrix6x13 C3D8I::enhanced_strain_matrix(const EnhancedModes& modes) {
    Matrix6x13 G = Matrix6x13::Zero();

    for (Index mode = 0; mode < n_modes; ++mode) {
        const Mat3& H = modes[mode];

        G(0, mode) = H(0, 0);
        G(1, mode) = H(1, 1);
        G(2, mode) = H(2, 2);
        G(3, mode) = H(1, 2) + H(2, 1);
        G(4, mode) = H(0, 2) + H(2, 0);
        G(5, mode) = H(0, 1) + H(1, 0);
    }

    return G;
}

/**
 * Builds the Green-Lagrange strain derivatives with respect to the enhanced
 * parameters.
 *
 * For one deformation-gradient variation H_alpha, the variation of
 *
 *     E = 1/2 (F^T F - I)
 *
 * is sym(F^T H_alpha). The resulting engineering-shear components are stored in
 * the same six-component ordering used by the volume material interface.
 *
 * @param deformation_gradient Current enhanced deformation gradient F.
 * @param modes Deformation-gradient variations dF/dalpha for all local modes.
 * @return Green-Lagrange B matrix with one column per enhanced parameter.
 */
C3D8I::Matrix6x13 C3D8I::enhanced_green_lagrange_matrix(
    const Mat3&          deformation_gradient,
    const EnhancedModes& modes
) {
    Matrix6x13 G = Matrix6x13::Zero();

    for (Index mode = 0; mode < n_modes; ++mode) {
        const Mat3 A = deformation_gradient.transpose() * modes[mode];

        G(0, mode) = A(0, 0);
        G(1, mode) = A(1, 1);
        G(2, mode) = A(2, 2);
        G(3, mode) = A(1, 2) + A(2, 1);
        G(4, mode) = A(0, 2) + A(2, 0);
        G(5, mode) = A(0, 1) + A(1, 0);
    }

    return G;
}

/**
 * Assembles the coupled linear nodal/enhanced stiffness blocks.
 *
 * The reference-configuration integration forms
 *
 *     Kuu = integral B^T C B dV0,
 *     Kua = integral B^T C G dV0,
 *     Kau = integral G^T C B dV0,
 *     Kaa = integral G^T C G dV0.
 *
 * The material is evaluated at zero linearized strain using committed state
 * only. No material trial history is modified. The returned local block system
 * is subsequently used by stiffness condensation, enhanced-state recovery and
 * thermal loading.
 *
 * @return Linear nodal/enhanced stiffness blocks before static condensation.
 */
C3D8I::EnhancedSystem C3D8I::assemble_linear_system() {
    EnhancedSystem system;

    // Collect the reference geometry and full C3D8 stiffness quadrature
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const auto&              scheme           = this->integration_scheme_stiffness();
    const VolumeStrainLinearized zero_strain;

    // Integrate compatible, coupling and enhanced stiffness blocks over dV0
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto B = this->strain_displacement(dN_dX);
        const auto modes = enhanced_gradient_modes(
            reference_coords, point.r, point.s, point.t);
        const Matrix6x13 G = enhanced_strain_matrix(modes);

        // Evaluate the zero-strain constitutive tangent without advancing history
        const Index state_row = this->mp_index(ip);
        const Precision* old_state =
            &(*this->_model_data->material_state_old)(state_row, 0);

        VolumeStressCauchy zero_stress;
        Mat6 C;
        evaluate_material(
            point.r, point.s, point.t,
            zero_strain, old_state, nullptr,
            zero_stress, C);

        // Assemble the four blocks of the local linear EAS system
        const Precision measure = point.w * det0;

        system.kuu.noalias() += measure * B.transpose() * C * B;
        system.kua.noalias() += measure * B.transpose() * C * G;
        system.kau.noalias() += measure * G.transpose() * C * B;
        system.kaa.noalias() += measure * G.transpose() * C * G;
    }

    return system;
}

/**
 * Solves the stationary enhanced parameters for a linearized displacement
 * state.
 *
 * The compatible displacement produces the local enhanced residual
 *
 *     r_alpha = Kau u.
 *
 * Optional isotropic thermal free strain contributes the corresponding
 * negative enhanced thermal force. The 13 x 13 enhanced stiffness is then
 * solved locally, so the returned parameters never enter the global system.
 *
 * @param displacement Element translational displacement vector.
 * @param thermal_free_strain Optional scalar free strain at the eight nodes.
 * @return Stationary local enhanced parameters alpha.
 */
C3D8I::Vector13 C3D8I::solve_linear_modes(
    const Vector24&        displacement,
    const StaticVector<N>* thermal_free_strain
) {
    const EnhancedSystem system = assemble_linear_system();
    Vector13 residual = system.kau * displacement;

    // Prescribed isotropic free strain also acts on the local EAS equations.
    // Condensing it here is required for stress recovery and the equivalent
    // thermal nodal load to remain consistent with the C3D8I stiffness.
    if (thermal_free_strain != nullptr) {
        const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
        const auto& scheme = this->integration_scheme_stiffness();
        const VolumeStrainLinearized zero_strain;

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);

            Precision det0 = Precision(0);
            const auto dN_dX = this->shape_derivatives_reference(
                reference_coords, point.r, point.s, point.t, det0);
            (void) dN_dX;

            const EnhancedModes modes = enhanced_gradient_modes(
                reference_coords, point.r, point.s, point.t);
            const Matrix6x13 G = enhanced_strain_matrix(modes);

            const Precision free_value =
                this->shape_function(point.r, point.s, point.t).dot(*thermal_free_strain);
            const Vec6 free_strain {
                free_value, free_value, free_value,
                Precision(0), Precision(0), Precision(0)
            };

            const Index state_row = this->mp_index(ip);
            const Precision* old_state =
                &(*this->_model_data->material_state_old)(state_row, 0);

            VolumeStressCauchy zero_stress;
            Mat6 C;
            evaluate_material(
                point.r, point.s, point.t,
                zero_strain, old_state, nullptr,
                zero_stress, C);

            residual.noalias() -=
                point.w * det0 * G.transpose() * C * free_strain;
        }
    }

    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular incompatible-mode stiffness in element ", elem_id);

    return -solver.solve(residual);
}

/**
 * Assembles and statically condenses the linear C3D8I stiffness.
 *
 * The four local EAS blocks are reduced by
 *
 *     K = Kuu - Kua Kaa^-1 Kau.
 *
 * Only the resulting 24 x 24 nodal matrix is exposed to the global assembler.
 * The analytically symmetric result is symmetrized only to remove numerical
 * round-off.
 *
 * @param buffer Caller-owned contiguous storage for the element matrix.
 * @return Map onto the condensed nodal stiffness matrix.
 */
MapMatrix C3D8I::stiffness(Precision* buffer) {
    // Assemble the coupled local system and factor the enhanced block
    const EnhancedSystem system = assemble_linear_system();

    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular incompatible-mode stiffness in element ", elem_id);

    // Eliminate all thirteen enhanced parameters before global assembly
    Matrix24 condensed = system.kuu - system.kua * solver.solve(system.kau);

    // Linear elasticity is analytically symmetric; remove round-off asymmetry in
    // the same way as the common fully integrated solid implementation.
    condensed = (Precision(0.5) * (condensed + condensed.transpose())).eval();

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = condensed;
    return mapped;
}

/**
 * Computes compliance sensitivity to the three additional material-orientation
 * angles using the stationary C3D8I strain field.
 *
 * For fixed nodal displacement, the condensed element energy is stationary with
 * respect to the enhanced parameters. The envelope theorem therefore removes
 * all explicit d(alpha)/d(theta) terms and the orientation derivative can be
 * evaluated directly from
 *
 *     epsilon = B u + G alpha
 *
 * as
 *
 *     dJ/dtheta_i = integral epsilon^T (dC/dtheta_i) epsilon dV0.
 *
 * This is the same orientation derivative used by the common solid element, but
 * evaluated with the complete stationary C3D8I strain instead of the compatible
 * C3D8 strain. Constitutive history remains state-neutral.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element-domain field receiving the three angle derivatives.
 */
void C3D8I::compute_compliance_angle_derivative(
    Field& displacement,
    Field& result
) {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return;
    }

    // Build the additional material rotation and its three angle derivatives
    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name, "': material orientation requires 3 components");

    const Index row    = static_cast<Index>(this->elem_id);
    const Vec3  angles = angles_field->row_vec3(row);

    const Mat3 additional_rotation = cos::RectangularSystem::euler(
        angles(0),
        angles(1),
        angles(2)
    ).get_axes(Vec3::Zero());

    const std::array<Mat3, 3> additional_rotation_derivatives {
        cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2))
    };

    // Recover the stationary enhanced strain state for the current orientation
    const Vector24 u     = local_displacement(displacement);
    const Vector13 alpha = solve_linear_modes(u);

    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const auto&              scheme           = this->integration_scheme_stiffness();

    const Precision scaling    = this->element_stiffness_scale();
    Vec3            derivative = Vec3::Zero();

    // Integrate the orientation sensitivity of the complete stationary strain
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );
        const auto B = this->strain_displacement(dN_dX);
        const EnhancedModes modes = enhanced_gradient_modes(
            reference_coords,
            point.r,
            point.s,
            point.t
        );
        const Matrix6x13 G = enhanced_strain_matrix(modes);
        const Vec6 strain  = B * u + G * alpha;

        // Differentiate only the constitutive orientation transformation while
        // reading the committed material-point history.
        const Vec3 position_reference =
            this->interpolate<D>(reference_coords, point.r, point.s, point.t);
        const Index      state_row = this->mp_index(ip);
        const Precision* old_state =
            &(*this->_model_data->material_state_old)(state_row, 0);

        const auto tangent_derivatives = this->get_section()->tangent_rotation_derivatives(
            position_reference,
            additional_rotation,
            additional_rotation_derivatives,
            old_state,
            nullptr
        );

        for (Index angle = 0; angle < 3; ++angle) {
            derivative(angle) +=
                scaling
                * point.w
                * det0
                * strain.dot(tangent_derivatives[angle] * strain);
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

MapMatrix C3D8I::stiffness_geom(Precision* buffer, const Field& displacement) {
    return stiffness_geom(buffer, displacement, nullptr);
}

/**
 * Builds the condensed initial-stress operator for linear buckling.
 *
 * The supplied preload displacement first defines the stationary linear EAS
 * state
 *
 *     epsilon = B u + G alpha.
 *
 * Cauchy prestress from this state produces the four geometric block matrices
 *
 *     Guu, Gua, Gau, Gaa.
 *
 * Linear buckling uses the first-order prestress perturbation of the already
 * condensed elastic element. Differentiating the Schur complement
 *
 *     Kc(lambda)
 *       = Kuu + lambda Guu
 *       - (Kua + lambda Gua)
 *         (Kaa + lambda Gaa)^-1
 *         (Kau + lambda Gau)
 *
 * at lambda = 0 gives
 *
 *     Gc = Guu
 *        - Gua Kaa^-1 Kau
 *        - Kua Kaa^-1 Gau
 *        + Kua Kaa^-1 Gaa Kaa^-1 Kau.
 *
 * This keeps the global generalized buckling problem linear in the load factor
 * while accounting for the perturbation of all element-local enhanced
 * parameters. Constitutive history remains state-neutral.
 *
 * @param buffer Caller-owned contiguous storage for the geometric matrix.
 * @param displacement Global nodal displacement field defining prestress.
 * @param thermal_free_strain Optional element-nodal scalar thermal free strain.
 * @return Map onto the condensed 24 x 24 initial-stress matrix.
 */
MapMatrix C3D8I::stiffness_geom(
    Precision*   buffer,
    const Field& displacement,
    const Field* thermal_free_strain
) {
    // Validate optional thermal prestress data
    if (thermal_free_strain != nullptr) {
        logging::error(thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL
                       && thermal_free_strain->components == 1,
            "C3D8I: thermal free strain must be scalar ELEMENT_NODAL data");
    }

    // Collect reference geometry, displacement and optional nodal free strain
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const Vector24           u                = local_displacement(displacement);

    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain != nullptr) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) =
                (*thermal_free_strain)(
                    static_cast<Index>(this->elem_nodal_offset) + node,
                    0
                );
        }
    }

    // Recover the stationary enhanced preload state and elastic condensation blocks
    const Vector13 alpha = solve_linear_modes(
        u,
        thermal_free_strain != nullptr ? &nodal_thermal_strain : nullptr
    );
    const EnhancedSystem elastic = assemble_linear_system();

    Eigen::FullPivLU<Matrix13> solver(elastic.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular incompatible-mode stiffness in element ", elem_id);

    // Initial-stress blocks before condensation
    Matrix24    guu = Matrix24::Zero();
    Matrix24x13 gua = Matrix24x13::Zero();
    Matrix13x24 gau = Matrix13x24::Zero();
    Matrix13    gaa = Matrix13::Zero();

    const auto& scheme = this->integration_scheme_stiffness();

    // Reconstruct prestress and integrate all nodal/enhanced geometric blocks
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );
        const auto B = this->strain_displacement(dN_dX);
        const EnhancedModes modes = enhanced_gradient_modes(
            reference_coords,
            point.r,
            point.s,
            point.t
        );
        const Matrix6x13 G = enhanced_strain_matrix(modes);

        Vec6 strain_values = B * u + G * alpha;
        if (thermal_free_strain != nullptr) {
            const Precision free_value =
                this->shape_function(point.r, point.s, point.t).dot(nodal_thermal_strain);
            strain_values(0) -= free_value;
            strain_values(1) -= free_value;
            strain_values(2) -= free_value;
        }

        // Evaluate the state-neutral Cauchy prestress at the current material point
        const VolumeStrainLinearized element_strain(strain_values);
        const Index      state_row = this->mp_index(ip);
        const Precision* old_state =
            &(*this->_model_data->material_state_old)(state_row, 0);

        VolumeStressCauchy stress;
        Mat6               material_tangent;
        evaluate_material(
            point.r,
            point.s,
            point.t,
            element_strain,
            old_state,
            nullptr,
            stress,
            material_tangent
        );

        const Mat3      sigma   = stress.tensor();
        const Precision measure = point.w * det0;

        // Nodal-nodal geometric block
        for (Index node_a = 0; node_a < N; ++node_a) {
            const Vec3 dNa = dN_dX.row(node_a).transpose();

            for (Index node_b = 0; node_b < N; ++node_b) {
                const Vec3 dNb = dN_dX.row(node_b).transpose();
                const Precision coefficient = dNa.dot(sigma * dNb) * measure;

                for (Dim component = 0; component < D; ++component) {
                    guu(D * node_a + component,
                        D * node_b + component) += coefficient;
                }
            }

            // Nodal-enhanced and enhanced-nodal geometric coupling
            for (Dim component = 0; component < D; ++component) {
                for (Index mode = 0; mode < n_modes; ++mode) {
                    const Vec3 h = modes[mode].row(component).transpose();
                    const Precision coefficient = dNa.dot(sigma * h) * measure;

                    gua(D * node_a + component, mode) += coefficient;
                    gau(mode, D * node_a + component) += coefficient;
                }
            }
        }

        // Enhanced-enhanced geometric block
        for (Index mode_a = 0; mode_a < n_modes; ++mode_a) {
            for (Index mode_b = 0; mode_b < n_modes; ++mode_b) {
                Precision coefficient = Precision(0);

                for (Dim component = 0; component < D; ++component) {
                    const Vec3 ha = modes[mode_a].row(component).transpose();
                    const Vec3 hb = modes[mode_b].row(component).transpose();
                    coefficient += ha.dot(sigma * hb);
                }

                gaa(mode_a, mode_b) += measure * coefficient;
            }
        }
    }

    // Differentiate the elastic Schur complement with respect to prestress.
    // The first solve is the enhanced response to an arbitrary nodal perturbation.
    const Matrix13x24 alpha_u = solver.solve(elastic.kau);

    Matrix24 geometric =
        guu
        - gua * alpha_u
        - elastic.kua * solver.solve(gau)
        + elastic.kua * solver.solve(gaa * alpha_u);

    // The exact operator is symmetric for a conservative material formulation
    geometric = (Precision(0.5) * (geometric + geometric.transpose())).eval();

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = geometric;
    return mapped;
}

/**
 * Assembles the condensed equivalent nodal load from isotropic thermal strain.
 *
 * The compatible thermal force f_u and local enhanced force f_a are integrated
 * together and the latter is removed through the same Schur complement as the
 * mechanical stiffness:
 *
 *     f_th = f_u - K_ua K_aa^-1 f_a.
 *
 * Material history remains state-neutral because thermal load construction uses
 * only committed constitutive state.
 *
 * @param node_loads Global nodal load accumulator.
 * @param node_temp Scalar nodal temperature field.
 * @param ref_temp Stress-free reference temperature.
 */
void C3D8I::apply_tload(
    Field&       node_loads,
    const Field& node_temp,
    Precision    ref_temp
) {
    // Validate all fields and material properties required by thermal loading
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
        "C3D8I: thermal loading requires a scalar nodal temperature field");
    logging::error(node_loads.domain == FieldDomain::NODE && node_loads.components >= D,
        "C3D8I: thermal loading requires nodal load storage with three components");
    logging::error(material() != nullptr && material()->has_elasticity(),
        "C3D8I: thermal loading requires elasticity in element ", elem_id);
    logging::error(material()->has_thermal_expansion(),
        "C3D8I: material has no thermal expansion in element ", elem_id);

    // Collect finite nodal temperatures and replace undefined values by T_ref
    StaticVector<N> temperatures = StaticVector<N>::Zero();
    for (Index node = 0; node < N; ++node) {
        const Precision value = node_temp(static_cast<Index>(node_ids[node]), 0);
        temperatures(node) = std::isfinite(value) ? value : ref_temp;
    }

    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();
    const VolumeStrainLinearized zero_strain;

    Vector24 compatible_load = Vector24::Zero();
    Vector13 enhanced_load   = Vector13::Zero();

    // Integrate compatible and enhanced thermal-force blocks over dV0
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto B = this->strain_displacement(dN_dX);
        const EnhancedModes modes = enhanced_gradient_modes(
            reference_coords, point.r, point.s, point.t);
        const Matrix6x13 G = enhanced_strain_matrix(modes);

        const Precision temperature =
            this->shape_function(point.r, point.s, point.t).dot(temperatures);
        const Precision free_value =
            material()->get_thermal_expansion() * (temperature - ref_temp);
        const Vec6 free_strain {
            free_value, free_value, free_value,
            Precision(0), Precision(0), Precision(0)
        };

        const Index state_row = this->mp_index(ip);
        const Precision* old_state =
            &(*this->_model_data->material_state_old)(state_row, 0);

        VolumeStressCauchy zero_stress;
        Mat6 C;
        evaluate_material(
            point.r, point.s, point.t,
            zero_strain, old_state, nullptr,
            zero_stress, C);

        const Vec6 thermal_stress = C * free_strain;
        const Precision measure = point.w * det0;
        compatible_load.noalias() += measure * B.transpose() * thermal_stress;
        enhanced_load.noalias() += measure * G.transpose() * thermal_stress;
    }

    // Condense the enhanced thermal force with the same Kua/Kaa blocks as K
    const EnhancedSystem system = assemble_linear_system();
    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular incompatible-mode stiffness in element ", elem_id);

    const Vector24 condensed_load =
        compatible_load - system.kua * solver.solve(enhanced_load);
    assemble_local_force(node_loads, condensed_load);
}

/**
 * Prepares the fixed geometry used throughout one local enhanced-state solve.
 *
 * Reference derivatives and enhanced modes are evaluated once at each full
 * C3D8 integration point. The compatible deformation gradient is formed as
 * F_c = x^T dN/dX; it remains fixed while the local enhanced parameters change.
 * This evaluation-local data is also reused by final assembly or recovery and
 * carries no constitutive history across independent trial configurations.
 *
 * @param reference_coords Global nodal positions in the reference configuration.
 * @param current_coords Global nodal positions in the supplied current configuration.
 * @return Eight pointwise geometries with complete reference-volume weights.
 */
C3D8I::NonlinearPoints C3D8I::nonlinear_points(
    const StaticMatrix<N, D>& reference_coords,
    const StaticMatrix<N, D>& current_coords
) {
    NonlinearPoints points;
    const auto& scheme = this->integration_scheme_stiffness();

    // Preserve the full eight-point rule inherited from C3D8
    logging::error(scheme.count() == points.size(),
        "C3D8I: nonlinear evaluation requires eight integration points");

    // Prepare immutable geometry before Newton iterates on the enhanced parameters
    for (Index ip = 0; ip < points.size(); ++ip) {
        const auto quadrature = scheme.get_point(ip);
        auto& point = points[ip];

        Precision det0 = Precision(0);
        point.natural     = Vec3(quadrature.r, quadrature.s, quadrature.t);
        point.derivatives = this->shape_derivatives_reference(
            reference_coords, quadrature.r, quadrature.s, quadrature.t, det0);
        point.compatible  = current_coords.transpose() * point.derivatives;
        point.measure     = quadrature.w * det0;
        point.modes       = enhanced_gradient_modes(
            reference_coords, quadrature.r, quadrature.s, quadrature.t);

        // Reject an inverted compatible configuration before any local material query
        const Precision determinant = point.compatible.determinant();
        logging::error(std::isfinite(determinant) && determinant > Precision(0),
            "C3D8I: non-positive compatible deformation gradient in element ", elem_id,
            "\ndet(F_c): ", determinant);
    }
    return points;
}

/**
 * Assembles the coupled nodal/enhanced residual and consistent finite-strain
 * tangent before static condensation.
 *
 * The enhanced tensors are defined in the reference configuration and enter the
 * deformation gradient through the right multiplicative map
 *
 *     M     = I + sum_m alpha_m H_m,
 *     F_bar = F_c M.
 *
 * A superposed spatial rigid rotation therefore gives
 *
 *     F_c -> Q F_c,
 *     F_bar -> Q F_bar,
 *
 * without changing the admissible enhanced mode space. This keeps the local
 * stationarity problem and the condensed element response objective.
 *
 * For fixed enhanced parameters, the nodal deformation-gradient variation is
 *
 *     delta F_u = delta F_c M.
 *
 * For one enhanced parameter it is
 *
 *     delta F_alpha = F_c H_alpha.
 *
 * The corresponding Green-Lagrange B matrices provide the material tangent
 * blocks. The stress-dependent geometric blocks use the same first variations.
 * Mixed nodal/enhanced derivatives additionally contain
 *
 *     d2 F / (du dalpha) = dF_c/du H_alpha,
 *
 * which contributes the second mixed geometric term required by the consistent
 * tangent.
 *
 * Material history is read from the committed state at every integration point.
 * A writable trial-state row is supplied only for the final physical nonlinear
 * evaluation after the local enhanced Newton solve has converged.
 *
 * @param points Fixed reference/current geometry at the eight integration points.
 * @param alpha Current element-local enhanced parameters.
 * @param write_material_state Write constitutive trial state for this evaluation.
 * @param assemble_global_blocks Also assemble nodal residual and nodal coupling blocks.
 * @param assemble_tangent Construct local/global tangents. False selects the
 *                         final nodal residual-only path after local convergence.
 * @return Coupled nodal/enhanced residual and tangent blocks before condensation.
 */
C3D8I::EnhancedSystem C3D8I::assemble_nonlinear_system(
    const NonlinearPoints& points,
    const Vector13&        alpha,
    bool                   write_material_state,
    bool                   assemble_global_blocks,
    bool                   assemble_tangent
) {
    EnhancedSystem system;

    // Integrate all residual and tangent blocks over the reference element volume
    for (Index ip = 0; ip < points.size(); ++ip) {
        const auto& point        = points[ip];
        const auto& dN_dX        = point.derivatives;
        const auto& compatible_F = point.compatible;
        const auto& modes        = point.modes;

        // Build the element-local right multiplicative enhancement
        Mat3 enhancement = Mat3::Identity();
        for (Index mode = 0; mode < n_modes; ++mode) {
            enhancement.noalias() += alpha(mode) * modes[mode];
        }

        const Mat3 F = compatible_F * enhancement;
        const Precision determinant = F.determinant();
        logging::error(std::isfinite(determinant) && determinant > Precision(0),
            "C3D8I: non-positive enhanced deformation gradient in element ", elem_id,
            "\ndet(F): ", determinant);

        // Transform the compatible nodal deformation-gradient variations through
        // the same right enhancement used by the total deformation gradient
        const StaticMatrix<N, D> enhanced_shape_derivatives = dN_dX * enhancement;

        // Evaluate the work-conjugate PK2 material response in the reference configuration
        const VolumeStrainGreenLagrange strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        const Index      state_row = this->mp_index(ip);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        Precision*       new_state = write_material_state
            ? &(*this->_model_data->material_state_new)(state_row, 0)
            : nullptr;

        VolumeStressPK2 stress;
        Mat6            C;
        evaluate_material(
            point.natural(0),
            point.natural(1),
            point.natural(2),
            strain,
            old_state,
            new_state,
            stress,
            assemble_tangent ? &C : nullptr
        );

        const Precision measure      = point.measure;
        const Vec6      stress_voigt = stress.voigt();

        // Final force-only evaluations need neither enhanced nor nodal Hessians
        if (!assemble_tangent) {
            const auto Bu = this->green_lagrange_strain_displacement(enhanced_shape_derivatives, F);
            system.ru.noalias() += measure * Bu.transpose() * stress_voigt;
            continue;
        }

        // Enhanced variations and their stress contractions are reused across all
        // enhanced-enhanced and nodal-enhanced geometric tangent entries.
        const Mat3 S = stress.tensor();
        EnhancedModes alpha_variations;
        EnhancedModes alpha_stress_variations;
        for (Index mode = 0; mode < n_modes; ++mode) {
            alpha_variations[mode]        = compatible_F * modes[mode];
            alpha_stress_variations[mode] = alpha_variations[mode] * S;
        }
        const Matrix6x13 Ba = enhanced_green_lagrange_matrix(F, alpha_variations);
        const Matrix6x13 CBa = measure * C * Ba;

        // Assemble the enhanced residual and its material tangent. These blocks
        // are required by every local alpha Newton iteration.
        system.ra.noalias()  += measure * Ba.transpose() * stress_voigt;
        system.kaa.noalias() += Ba.transpose() * CBa;

        // The local Newton solve does not need any nodal tangent or residual block
        if (assemble_global_blocks) {
            const auto Bu = this->green_lagrange_strain_displacement(enhanced_shape_derivatives, F);
            const StaticMatrix<6, ndof> CBu = measure * C * Bu;
            system.ru.noalias()  += measure * Bu.transpose() * stress_voigt;
            system.kuu.noalias() += Bu.transpose() * CBu;
            system.kua.noalias() += Bu.transpose() * CBa;
            system.kau.noalias() += Ba.transpose() * CBu;
        }

        // Add the enhanced-enhanced geometric tangent from
        // delta F_alpha = F_c H_alpha
        for (Index mode_a = 0; mode_a < n_modes; ++mode_a) {
            for (Index mode_b = 0; mode_b < n_modes; ++mode_b) {
                const Precision coefficient =
                    alpha_variations[mode_a].cwiseProduct(alpha_stress_variations[mode_b]).sum();
                system.kaa(mode_a, mode_b) += measure * coefficient;
            }
        }

        if (!assemble_global_blocks) {
            continue;
        }

        // Add nodal-nodal and nodal-enhanced geometric tangent contributions.
        // The mixed block includes the second derivative of F_c M with respect
        // to one nodal and one enhanced parameter.
        const Mat3 FS = F * S;

        for (Index node_a = 0; node_a < N; ++node_a) {
            const Vec3 dNa          = dN_dX.row(node_a).transpose();
            const Vec3 enhanced_dNa = enhanced_shape_derivatives.row(node_a).transpose();

            // Nodal-nodal geometric tangent from delta F_u = delta F_c M
            for (Index node_b = 0; node_b < N; ++node_b) {
                const Vec3 enhanced_dNb =
                    enhanced_shape_derivatives.row(node_b).transpose();
                const Precision coefficient =
                    enhanced_dNa.dot(S * enhanced_dNb) * measure;

                for (Dim component = 0; component < D; ++component) {
                    system.kuu(D * node_a + component,
                               D * node_b + component) += coefficient;
                }
            }

            // Nodal-enhanced geometric tangent including d2F/(du dalpha)
            for (Index mode = 0; mode < n_modes; ++mode) {
                const Vec3 first_variation = alpha_stress_variations[mode] * enhanced_dNa;
                const Vec3 mixed_variation = FS * (modes[mode].transpose() * dNa);
                const Vec3 coefficient     = measure * (first_variation + mixed_variation);

                for (Dim component = 0; component < D; ++component) {
                    system.kua(D * node_a + component, mode) += coefficient(component);
                    system.kau(mode, D * node_a + component) += coefficient(component);
                }
            }
        }
    }

    return system;
}

/**
 * Solves the thirteen finite-strain enhanced stationarity equations.
 *
 * The nodal reference and current coordinates remain fixed while Newton
 * iterates only on the local enhanced parameters. Every iteration assembles
 * r_alpha and K_alpha_alpha from committed constitutive history without writing
 * trial state. The converged alpha is subsequently used by the physical
 * nonlinear tangent evaluation, which performs the actual material-state update.
 *
 * @param points Fixed reference/current geometry at the eight integration points.
 * @return Converged element-local enhanced parameters.
 */
C3D8I::Vector13 C3D8I::solve_nonlinear_modes(const NonlinearPoints& points) {
    // Initialize the element-local enhanced state for the current trial geometry
    Vector13 alpha = Vector13::Zero();

    constexpr Index     max_iterations = 20;
    constexpr Precision tolerance      = Precision(1e-10);

    // Enforce stationarity of the thirteen local enhanced equations
    for (Index iteration = 0; iteration < max_iterations; ++iteration) {
        const EnhancedSystem system = assemble_nonlinear_system(
            points, alpha, false, false, true);

        Eigen::FullPivLU<Matrix13> solver(system.kaa);
        logging::error(solver.isInvertible(),
            "C3D8I: singular local EAS tangent in element ", elem_id);

        // Apply the local Newton correction without exposing alpha globally
        const Vector13 delta = solver.solve(-system.ra);
        alpha += delta;

        // Scale the stopping criterion with the current enhanced-parameter level
        const Precision scale = Precision(1) + alpha.cwiseAbs().maxCoeff();
        if (delta.cwiseAbs().maxCoeff() <= tolerance * scale) {
            return alpha;
        }
    }

    logging::error(false,
        "C3D8I: local incompatible-mode Newton did not converge in element ", elem_id);
    return alpha;
}

C3D8I::Vector24 C3D8I::local_displacement(const Field& displacement) {
    const StaticMatrix<N, D> local = this->nodal_data<D>(displacement);
    Vector24 result = Vector24::Zero();

    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = local(node, dof);
        }
    }

    return result;
}

void C3D8I::assemble_local_force(Field& node_forces, const Vector24& local_force) {
    logging::error(node_forces.domain == FieldDomain::NODE,
        "C3D8I: internal force output must use NODE domain");
    logging::error(node_forces.components >= D,
        "C3D8I: internal force output requires at least three components");

    for (Index node = 0; node < N; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);
        for (Dim dof = 0; dof < D; ++dof) {
            node_forces(node_id, dof) += local_force(D * node + dof);
        }
    }
}

/**
 * Evaluates the nonlinear C3D8I internal force and consistent condensed tangent.
 *
 * The nodal trial configuration first defines the compatible deformation. A
 * local Newton solve enforces r_alpha = 0 for the thirteen enhanced parameters.
 * The converged state is then re-evaluated with writable constitutive trial
 * history, after which the local tangent is reduced by
 *
 *     K = Kuu - Kua Kaa^-1 Kau.
 *
 * The internal force is already stationary with respect to alpha and is
 * scattered directly into the global nodal accumulator.
 *
 * @param buffer Optional caller-owned storage for the 24 x 24 tangent. A null
 *               pointer requests internal force only.
 * @param nodal_forces Global nodal internal-force accumulator.
 * @param displacement Global nodal displacement field of the current trial state.
 * @return Map onto the condensed tangent, or an empty map when buffer is null.
 */
MapMatrix C3D8I::stiffness_tangent(
    Precision*   buffer,
    NodeData&    nodal_forces,
    const Field& displacement
) {
    // Construct the compatible current geometry from the global displacement field
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const StaticMatrix<N, D> local_u          = this->nodal_data<D>(displacement);
    const StaticMatrix<N, D> current_coords   = reference_coords + local_u;

    // Enforce stationarity of the element-local enhanced parameters
    const NonlinearPoints points = nonlinear_points(reference_coords, current_coords);
    const Vector13        alpha  = solve_nonlinear_modes(points);

    // Re-evaluate the converged state and write exactly this constitutive trial state
    const EnhancedSystem system = assemble_nonlinear_system(
        points, alpha, true, true, buffer != nullptr);

    // Scatter the stationary nodal internal force before optional tangent assembly
    assemble_local_force(nodal_forces, system.ru);

    if (buffer == nullptr) {
        return MapMatrix(nullptr, 0, 0);
    }

    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular converged local EAS tangent in element ", elem_id);

    // Condense the local enhanced unknowns from the global Newton tangent
    const Matrix24 condensed = system.kuu - system.kua * solver.solve(system.kau);

    MapMatrix mapped{buffer, ndof, ndof};
    mapped = condensed;
    return mapped;
}

void C3D8I::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    bool             use_green_lagrange_nl
) {
    compute_stress_strain(
        strain, stress, displacement, rst, offset,
        use_green_lagrange_nl, nullptr);
}

/**
 * Recovers physical strain and Cauchy stress from the stationary enhanced state.
 *
 * Linear recovery solves the condensed enhanced parameters from the supplied
 * displacement and optional thermal free strain. Nonlinear recovery solves the
 * objective multiplicative enhanced state
 *
 *     F_bar = F_c (I + sum_m alpha_m H_m)
 *
 * and evaluates Green-Lagrange strain with PK2 constitutive stress before
 * pushing the result forward to Cauchy stress.
 *
 * Constitutive history is never advanced during recovery. Integration-point
 * values are either written directly or extrapolated to the natural element
 * nodes using the inherited C3D8 extrapolation operator.
 *
 * @param strain Optional strain output field.
 * @param stress Optional Cauchy-stress output field.
 * @param displacement Global nodal displacement field.
 * @param rst Requested natural output coordinates.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select finite-strain or linearized recovery.
 * @param thermal_free_strain Optional scalar element-nodal thermal free strain.
 */
void C3D8I::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    bool             use_green_lagrange_nl,
    const Field*     thermal_free_strain
) {
    // Validate requested output and supported thermal/nonlinear combination
    logging::error(strain != nullptr || stress != nullptr,
        "C3D8I: stress/strain recovery requires at least one output field");
    logging::error(!use_green_lagrange_nl || thermal_free_strain == nullptr,
        "C3D8I: nonlinear thermal recovery is not supported");
    if (thermal_free_strain != nullptr) {
        logging::error(thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1,
            "C3D8I: thermal free strain must be scalar ELEMENT_NODAL data");
    }

    // Classify output coordinates as integration-point or nodal recovery
    const auto&     scheme = this->integration_scheme_stiffness();
    const RowMatrix ip_rst = this->stress_strain_ip_rst();
    const bool output_at_ip =
        rst.rows() == ip_rst.rows() && rst.leftCols(3).isApprox(ip_rst);
    const bool output_at_nodes = rst.rows() == static_cast<Eigen::Index>(N);

    logging::error(output_at_ip || output_at_nodes,
        "C3D8I: stress/strain output must use integration points or element nodes");

    // Collect compatible reference/current geometry and element displacement
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const StaticMatrix<N, D> local_u          = this->nodal_data<D>(displacement);
    const StaticMatrix<N, D> current_coords   = reference_coords + local_u;
    const Vector24           u                = local_displacement(displacement);

    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain != nullptr) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) =
                (*thermal_free_strain)(
                    static_cast<Index>(this->elem_nodal_offset) + node, 0);
        }
    }

    // Reconstruct the stationary enhanced state for the selected kinematics
    NonlinearPoints points;
    Vector13        alpha;
    if (use_green_lagrange_nl) {
        points = nonlinear_points(reference_coords, current_coords);
        alpha  = solve_nonlinear_modes(points);
    } else {
        alpha = solve_linear_modes(
            u,
            thermal_free_strain != nullptr ? &nodal_thermal_strain : nullptr
        );
    }

    RowMatrix ip_strain = RowMatrix::Zero(scheme.count(), 6);
    RowMatrix ip_stress = RowMatrix::Zero(scheme.count(), 6);

    // Evaluate state-neutral strain and stress at every constitutive point
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const Index state_row = this->mp_index(ip);
        const Precision* old_state =
            &(*this->_model_data->material_state_old)(state_row, 0);

        if (!use_green_lagrange_nl) {
            Precision det0 = Precision(0);
            const auto dN_dX = this->shape_derivatives_reference(
                reference_coords, point.r, point.s, point.t, det0);
            const EnhancedModes modes = enhanced_gradient_modes(
                reference_coords, point.r, point.s, point.t);
            const auto B = this->strain_displacement(dN_dX);
            const Matrix6x13 G = enhanced_strain_matrix(modes);
            const Vec6 strain_values = B * u + G * alpha;
            Vec6 mechanical_strain = strain_values;

            if (thermal_free_strain != nullptr) {
                const Precision free_value =
                    this->shape_function(point.r, point.s, point.t)
                        .dot(nodal_thermal_strain);
                mechanical_strain(0) -= free_value;
                mechanical_strain(1) -= free_value;
                mechanical_strain(2) -= free_value;
            }

            const VolumeStrainLinearized element_strain(mechanical_strain);
            VolumeStressCauchy element_stress;
            Mat6 material_tangent;
            evaluate_material(
                point.r, point.s, point.t,
                element_strain, old_state, nullptr,
                element_stress, material_tangent);

            ip_strain.row(ip) = strain_values.transpose();
            ip_stress.row(ip) = element_stress.voigt().transpose();
            continue;
        }

        // Reconstruct the same objective multiplicative enhanced state used by
        // nonlinear residual and tangent assembly
        const Mat3&          compatible_F = points[ip].compatible;
        const EnhancedModes& modes        = points[ip].modes;

        Mat3 enhancement = Mat3::Identity();
        for (Index mode = 0; mode < n_modes; ++mode) {
            enhancement.noalias() += alpha(mode) * modes[mode];
        }

        const Mat3 F = compatible_F * enhancement;
        logging::error(std::isfinite(F.determinant()) && F.determinant() > Precision(0),
            "C3D8I: non-positive enhanced deformation gradient during recovery in element ",
            elem_id);

        const VolumeStrainGreenLagrange element_strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        evaluate_material(
            point.r, point.s, point.t,
            element_strain, old_state, nullptr,
            second_pk, nullptr);
        const VolumeStressCauchy element_stress = second_pk.to_cauchy(F);

        ip_strain.row(ip) = element_strain.voigt().transpose();
        ip_stress.row(ip) = element_stress.voigt().transpose();
    }

    // Write integration-point output directly when requested
    if (output_at_ip) {
        for (Eigen::Index row = 0; row < rst.rows(); ++row) {
            const Index global_row = static_cast<Index>(offset + row);
            for (Dim component = 0; component < 6; ++component) {
                if (strain) (*strain)(global_row, component) = ip_strain(row, component);
                if (stress) (*stress)(global_row, component) = ip_stress(row, component);
            }
        }
        return;
    }

    // Extrapolate constitutive-point values to the eight natural element nodes
    const RowMatrix& E = this->extrapolation_matrix();
    const RowMatrix nodal_strain = E * ip_strain;
    const RowMatrix nodal_stress = E * ip_stress;

    for (Eigen::Index row = 0; row < rst.rows(); ++row) {
        const Index global_row = static_cast<Index>(offset + row);
        for (Dim component = 0; component < 6; ++component) {
            if (strain) (*strain)(global_row, component) = nodal_strain(row, component);
            if (stress) (*stress)(global_row, component) = nodal_stress(row, component);
        }
    }
}

} // namespace fem::model
