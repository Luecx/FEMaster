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
 * Every expansion point uses multiplicative reference enhancement and a local
 * Newton solve of the thirteen stationarity equations. Static condensation
 * supplies the complete nodal tangent. Linear mechanics uses this same operator
 * at zero nodal displacement; finite mechanics moves the expansion point.
 *
 * The finite-strain kinematics use
 *
 * F_bar = F_compatible (I + sum_m alpha_m H_m),
 *
 * so a superposed spatial rigid rotation acts as
 * F_bar -> Q F_bar without changing the local enhanced mode space.
 *
 * Thermal loading, stress recovery and perturbation geometric stiffness reconstruct
 * the same stationary enhanced base state. Material history remains state-neutral in
 * all auxiliary paths; only evaluate() with update_state=true writes trial state.
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
 * alpha_i * xi_i
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
        "C3D8I: invalid reference determinant in element ", elem_id, "\ndet(J): ", j);
    logging::error(std::isfinite(j0) && j0 > Precision(0),
        "C3D8I: invalid center reference determinant in element ", elem_id, "\ndet(J0): ", j0);

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
 * Builds the Green-Lagrange strain derivatives with respect to the enhanced
 * parameters.
 *
 * For one deformation-gradient variation H_alpha, the variation of
 *
 * E = 1/2 (F^T F - I)
 *
 * is sym(F^T H_alpha). The resulting engineering-shear components are stored in
 * the same six-component ordering used by the volume material interface.
 *
 * @param deformation_gradient Current enhanced deformation gradient F.
 * @param modes Deformation-gradient variations dF/dalpha for all local modes.
 * @return Green-Lagrange B matrix with one column per enhanced parameter.
 */
C3D8I::Matrix6x13 C3D8I::enhanced_green_lagrange_matrix(
    const Mat3& deformation_gradient,
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
 * Queries the complete EAS block tangent at zero nodal displacement.
 *
 * The reference operator is the finite formulation evaluated at its stationary
 * enhanced state. Auxiliary queries read committed history and never write it.
 *
 * @return Coupled reference residual and tangent blocks before condensation.
 */
C3D8I::EnhancedSystem C3D8I::assemble_reference_system() {
    const auto reference = this->node_coords_reference();
    const auto points    = nonlinear_points(reference, reference);
    const auto alpha     = solve_nonlinear_modes(points);
    return assemble_nonlinear_system(points, alpha, false, true, true);
}

/**
 * Solves the enhanced parameter increment for a reference affine displacement.
 *
 * The finite reference tangent gives delta_alpha = -Kaa^-1 Kau u.
 * This auxiliary query supplies the reference compliance sensitivity.
 *
 * @param displacement Requested element translation increment from zero.
 * @return Enhanced parameter increment about the stationary reference state.
 */
C3D8I::Vector13 C3D8I::solve_linear_modes(const Vector24& displacement) {
    const auto reference = this->node_coords_reference();
    const auto points    = nonlinear_points(reference, reference);
    const auto alpha     = solve_nonlinear_modes(points);
    const auto system    = assemble_nonlinear_system(points, alpha, false, true, true);

    // Differentiate the stationary local equation using its complete tangent
    const Vector13 residual = system.kau * displacement;
    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular reference enhanced tangent in element ", elem_id);
    return -solver.solve(residual);
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
 * epsilon = B u + G alpha
 *
 * as
 *
 * dJ/dtheta_i = integral epsilon^T (dC/dtheta_i) epsilon dV0.
 *
 * This is the same orientation derivative used by the common solid element, but
 * evaluated with the complete stationary C3D8I strain instead of the compatible
 * C3D8 strain. Constitutive history remains state-neutral.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element-domain field receiving the three angle derivatives.
 */
void C3D8I::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return;
    }

    // Build the additional material rotation and its three angle derivatives
    auto angles_field                       = this->_model_data->material_orientation;
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

        Precision det0            = Precision(0);
        const auto dN_dX = this->shape_derivatives_reference(reference_coords, point.r, point.s, point.t, det0);
        const auto B              = this->strain_displacement(dN_dX);
        const EnhancedModes modes = enhanced_gradient_modes(reference_coords, point.r, point.s, point.t);
        const Matrix6x13 G        = enhanced_green_lagrange_matrix(Mat3::Identity(), modes);
        const Vec6 strain         = B * u + G * alpha;

        // Differentiate only the constitutive orientation transformation while
        // reading the committed material-point history.
        const Vec3 position_reference = this->interpolate<D>(reference_coords, point.r, point.s, point.t);
        const Index      state_row    = this->mp_index(ip);
        const Precision* old_state    = &(*this->_model_data->material_state_old)(state_row, 0);

        const auto tangent_derivatives = this->get_section()->tangent_rotation_derivatives(
            position_reference,
            additional_rotation,
            additional_rotation_derivatives,
            old_state,
            nullptr
        );

        for (Index angle = 0; angle < 3; ++angle) {
            derivative(angle) += scaling
                * point.w
                * det0
                * strain.dot(tangent_derivatives[angle] * strain);
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

/**
 * Builds the condensed geometric stiffness generated by a perturbation from u0
 * to u.
 *
 * The supplied nonlinear points, enhanced parameters and block system describe
 * the stationary base state u0. The nodal perturbation
 *
 *     Delta u = u - u0
 *
 * first induces the compatible/enhanced strain increment. Local stationarity
 * gives the corresponding enhanced-parameter increment
 *
 *     Delta alpha = -Kaa^-1 (Kau Delta u - f_alpha,th).
 *
 * The constitutive tangent at u0 then produces the PK2 stress increment
 *
 *     Delta S = C0 Delta E - Delta S_th.
 *
 * Only this stress increment is contracted with the second kinematic
 * derivatives. The stress already present at u0 remains part of the complete
 * base tangent and is not repeated here.
 *
 * Before condensation the perturbation geometric blocks are
 *
 *     Guu, Gua, Gau, Gaa.
 *
 * Differentiating the complete base-state Schur complement
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
 * @param buffer Caller-owned storage for the condensed geometric matrix.
 * @param points Finite kinematics at the base state u0.
 * @param alpha Stationary enhanced parameters at u0.
 * @param system Complete coupled tangent blocks at u0.
 * @param displacement_increment Nodal perturbation Delta u = u - u0.
 * @param thermal_free_strain Optional reference thermal free-strain perturbation.
 * @return Map onto the condensed 24 x 24 perturbation geometric matrix.
 */
MapMatrix C3D8I::stiffness_geom(
    Precision*             buffer,
    const NonlinearPoints& points,
    const Vector13&        alpha,
    const EnhancedSystem&  system,
    const Vector24&        displacement_increment,
    const Field*           thermal_free_strain
) {
    // Gather the optional scalar thermal free strain in element-node ordering.
    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) = (*thermal_free_strain)(this->elem_nodal_offset + node, 0);
        }
    }

    // Differentiate the enhanced stationarity equation at the base state.
    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular enhanced base tangent in element ", elem_id);

    Vector13 enhanced_residual = system.kau * displacement_increment;
    if (thermal_free_strain) {
        enhanced_residual -= thermal_force(points, alpha, nodal_thermal_strain).second;
    }
    const Vector13 alpha_increment = -solver.solve(enhanced_residual);

    StaticMatrix<N, D> local_delta = StaticMatrix<N, D>::Zero();
    for (Index node = 0; node < N; ++node) {
        local_delta.row(node) = displacement_increment.template segment<D>(D * node).transpose();
    }

    Matrix24    guu = Matrix24::Zero();
    Matrix24x13 gua = Matrix24x13::Zero();
    Matrix13x24 gau = Matrix13x24::Zero();
    Matrix13    gaa = Matrix13::Zero();

    // Recover Delta S at every material point and contract it with the same
    // second kinematic derivatives used by the complete finite tangent.
    for (Index ip = 0; ip < points.size(); ++ip) {
        const auto& point = points[ip];

        Mat3 enhancement       = Mat3::Identity();
        Mat3 delta_enhancement = Mat3::Zero();
        for (Index mode = 0; mode < n_modes; ++mode) {
            enhancement       += alpha(mode) * point.modes[mode];
            delta_enhancement += alpha_increment(mode) * point.modes[mode];
        }

        const Mat3 F = point.compatible * enhancement;
        const StaticMatrix<N, D> enhanced_shape_derivatives = point.derivatives * enhancement;

        // Compatible and enhanced perturbations contribute to the same
        // deformation-gradient increment about u0.
        const Mat3 delta_compatible = local_delta.transpose() * point.derivatives;
        const Mat3 delta_F = delta_compatible * enhancement + point.compatible * delta_enhancement;
        const Mat3 delta_E = Precision(0.5) * (F.transpose() * delta_F + delta_F.transpose() * F);

        Vec6 mechanical_increment = VolumeStrainGreenLagrange(delta_E).voigt();
        if (thermal_free_strain) {
            const Precision free = this->shape_function(point.natural(0), point.natural(1), point.natural(2))
                .dot(nodal_thermal_strain);
            mechanical_increment.head<3>().array() -= free;
        }

        const VolumeStrainGreenLagrange strain = VolumeStrainGreenLagrange::from_deformation_gradient(F);
        const Precision* old_state = &(*this->_model_data->material_state_old)(this->mp_index(ip), 0);

        VolumeStressPK2 stress_base;
        Mat6            material_tangent;
        evaluate_material(
            point.natural(0), point.natural(1), point.natural(2), strain,
            old_state, nullptr, stress_base, &material_tangent);

        const Mat3 stress_increment =
            VolumeStressPK2(Vec6(material_tangent * mechanical_increment)).tensor();

        // Enhanced first variations at the base state. Replacing S0 by Delta S
        // in the geometric contractions gives the perturbation block matrices.
        EnhancedModes alpha_variations;
        EnhancedModes alpha_stress_variations;
        for (Index mode = 0; mode < n_modes; ++mode) {
            alpha_variations[mode]        = point.compatible * point.modes[mode];
            alpha_stress_variations[mode] = alpha_variations[mode] * stress_increment;
        }

        const Precision measure = point.measure;

        for (Index mode_a = 0; mode_a < n_modes; ++mode_a) {
            for (Index mode_b = 0; mode_b < n_modes; ++mode_b) {
                const Precision coefficient =
                    alpha_variations[mode_a].cwiseProduct(alpha_stress_variations[mode_b]).sum();
                gaa(mode_a, mode_b) += measure * coefficient;
            }
        }

        const Mat3 F_stress_increment = F * stress_increment;

        for (Index node_a = 0; node_a < N; ++node_a) {
            const Vec3 dNa          = point.derivatives.row(node_a).transpose();
            const Vec3 enhanced_dNa = enhanced_shape_derivatives.row(node_a).transpose();

            for (Index node_b = 0; node_b < N; ++node_b) {
                const Vec3 enhanced_dNb     = enhanced_shape_derivatives.row(node_b).transpose();
                const Precision coefficient = measure * enhanced_dNa.dot(stress_increment * enhanced_dNb);

                for (Dim component = 0; component < D; ++component) {
                    guu(D * node_a + component, D * node_b + component) += coefficient;
                }
            }

            for (Index mode = 0; mode < n_modes; ++mode) {
                const Vec3 first_variation = alpha_stress_variations[mode] * enhanced_dNa;
                const Vec3 mixed_variation = F_stress_increment * (point.modes[mode].transpose() * dNa);
                const Vec3 coefficient     = measure * (first_variation + mixed_variation);

                for (Dim component = 0; component < D; ++component) {
                    gua(D * node_a + component, mode) += coefficient(component);
                    gau(mode, D * node_a + component) += coefficient(component);
                }
            }
        }
    }

    // Differentiate the complete base-state Schur complement with respect to
    // the perturbation stress.
    const Matrix13x24 alpha_u = solver.solve(system.kau);

    Matrix24 geometric = guu
        - gua * alpha_u
        - system.kua * solver.solve(gau)
        + system.kua * solver.solve(gaa * alpha_u);

    geometric = (Precision(0.5) * (geometric + geometric.transpose())).eval();

    MapMatrix mapped(buffer, ndof, ndof);
    mapped = geometric;
    return mapped;
}

/**
 * Integrates reference thermal-force increments in the finite EAS basis.
 *
 * At the reference stationary state, -C epsilon_th is an affine stress source.
 * Its compatible and enhanced force increments share the finite strain operators
 * used by the block tangent. No material history is modified.
 *
 * @param points Geometry at the reference expansion point.
 * @param alpha Stationary enhanced parameters at that point.
 * @param free_strain Scalar isotropic free strain at the eight element nodes.
 * @return Compatible and enhanced positive thermal-force vectors before condensation.
 */
std::pair<C3D8I::Vector24, C3D8I::Vector13> C3D8I::thermal_force(
    const NonlinearPoints& points,
    const Vector13&        alpha,
    const StaticVector<N>& free_strain
) {
    Vector24 force_u = Vector24::Zero();
    Vector13 force_a = Vector13::Zero();

    // Integrate the affine thermal source using the reference finite tangent
    for (Index ip = 0; ip < points.size(); ++ip) {
        const auto& point = points[ip];
        Mat3 enhancement = Mat3::Identity();
        for (Index mode = 0; mode < n_modes; ++mode) enhancement += alpha(mode) * point.modes[mode];
        const Mat3 F = point.compatible * enhancement;
        const auto strain = VolumeStrainGreenLagrange::from_deformation_gradient(F);
        const auto Bu = this->green_lagrange_strain_displacement(
            StaticMatrix<N, D>(point.derivatives * enhancement), F);
        EnhancedModes variations;
        for (Index mode = 0; mode < n_modes; ++mode) variations[mode] = point.compatible * point.modes[mode];
        const Matrix6x13 Ba = enhanced_green_lagrange_matrix(F, variations);

        const Precision* old_state = &(*this->_model_data->material_state_old)(this->mp_index(ip), 0);
        VolumeStressPK2 stress;
        Mat6            tangent;
        evaluate_material(
            point.natural(0), point.natural(1), point.natural(2), strain, old_state, nullptr, stress, &tangent);
        const Precision free = this->shape_function(point.natural(0), point.natural(1), point.natural(2)).dot(free_strain);
        const Vec6 thermal_stress = tangent * Vec6(free, free, free, 0, 0, 0);
        force_u.noalias() += point.measure * Bu.transpose() * thermal_stress;
        force_a.noalias() += point.measure * Ba.transpose() * thermal_stress;
    }
    return {force_u, force_a};
}

/**
 * Assembles reference thermal loading with the finite EAS Schur complement.
 *
 * Undefined temperatures use the stress-free reference temperature. The positive
 * thermal force is f_u - Kua Kaa^-1 f_a, matching the affine mechanical residual.
 *
 * @param node_loads Global nodal load accumulator.
 * @param node_temp Scalar nodal temperature field.
 * @param ref_temp Stress-free reference temperature.
 */
void C3D8I::apply_thermal_expansion_load(Field& node_loads, const Field& node_temp) {
    // Validate temperatures and form scalar nodal free strain
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
        "C3D8I: thermal loading requires a scalar nodal temperature field");
    auto mat = material();
    logging::error(mat != nullptr,
        "C3D8I: no material assigned to element ", elem_id);
    if (!mat->has_thermal_expansion()) {
        return;
    }

    const Precision ref_temp = mat->get_thermal_zero_temperature();
    const Precision alpha    = mat->get_thermal_expansion();

    StaticVector<N> free = StaticVector<N>::Zero();
    for (Index node = 0; node < N; ++node) {
        const Precision value = node_temp(node_ids[node], 0);
        const Precision temperature = std::isfinite(value) ? value : ref_temp;
        free(node) = alpha * (temperature - ref_temp);
    }

    // Condense the thermal source with the complete stationary reference tangent
    const auto reference = this->node_coords_reference();
    const auto points    = nonlinear_points(reference, reference);
    const auto alpha     = solve_nonlinear_modes(points);
    const auto system    = assemble_nonlinear_system(points, alpha, false, true, true);
    const auto thermal   = thermal_force(points, alpha, free);
    Eigen::FullPivLU<Matrix13> solver(system.kaa);
    logging::error(solver.isInvertible(),
        "C3D8I: singular reference enhanced tangent in element ", elem_id);
    assemble_local_force(node_loads, Vector24(thermal.first - system.kua * solver.solve(thermal.second)));
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
        auto& point           = points[ip];

        Precision det0    = Precision(0);
        point.natural     = Vec3(quadrature.r, quadrature.s, quadrature.t);
        point.derivatives = this->shape_derivatives_reference(
            reference_coords, quadrature.r, quadrature.s, quadrature.t, det0);
        point.compatible = current_coords.transpose() * point.derivatives;
        point.measure    = quadrature.w * det0;
        point.modes      = enhanced_gradient_modes(reference_coords, quadrature.r, quadrature.s, quadrature.t);

        // Reject an inverted compatible configuration before any local material query
        const Precision determinant = point.compatible.determinant();
        logging::error(std::isfinite(determinant) && determinant > Precision(0),
            "C3D8I: non-positive compatible deformation gradient in element ", elem_id, "\ndet(F_c): ", determinant);
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
 * M     = I + sum_m alpha_m H_m,
 * F_bar = F_c M.
 *
 * A superposed spatial rigid rotation therefore gives
 *
 * F_c -> Q F_c,
 * F_bar -> Q F_bar,
 *
 * without changing the admissible enhanced mode space. This keeps the local
 * stationarity problem and the condensed element response objective.
 *
 * For fixed enhanced parameters, the nodal deformation-gradient variation is
 *
 * delta F_u = delta F_c M.
 *
 * For one enhanced parameter it is
 *
 * delta F_alpha = F_c H_alpha.
 *
 * The corresponding Green-Lagrange B matrices provide the material tangent
 * blocks. The stress-dependent geometric blocks use the same first variations.
 * Mixed nodal/enhanced derivatives additionally contain
 *
 * d2 F / (du dalpha) = dF_c/du H_alpha,
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
 * final nodal residual-only path after local convergence.
 * @return Coupled nodal/enhanced residual and tangent blocks before condensation.
 */
C3D8I::EnhancedSystem C3D8I::assemble_nonlinear_system(
    const NonlinearPoints& points,
    const Vector13&        alpha,
    bool                   write_material_state,
    bool                   assemble_global_blocks,
    bool                   assemble_tangent,
    bool                   include_geometric,
    const StaticVector<N>* thermal_free_strain
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

        const Mat3 F                = compatible_F * enhancement;
        const Precision determinant = F.determinant();
        logging::error(std::isfinite(determinant) && determinant > Precision(0),
            "C3D8I: non-positive enhanced deformation gradient in element ", elem_id, "\ndet(F): ", determinant);

        // Transform the compatible nodal deformation-gradient variations through
        // the same right enhancement used by the total deformation gradient
        const StaticMatrix<N, D> enhanced_shape_derivatives = dN_dX * enhancement;

        // Evaluate the work-conjugate PK2 material response in the reference configuration
        const VolumeStrainGreenLagrange strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        Vec6 constitutive_strain_values = strain.voigt();
        if (thermal_free_strain) {
            const Precision free =
                this->shape_function(
                    point.natural(0), point.natural(1), point.natural(2))
                    .dot(*thermal_free_strain);
            constitutive_strain_values.head<3>().array() -= free;
        }
        const VolumeStrainGreenLagrange constitutive_strain(
            constitutive_strain_values);

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
            constitutive_strain,
            old_state,
            new_state,
            stress,
            assemble_tangent ? &C : nullptr
        );

        const Precision measure      = point.measure;
        const Vec6      stress_voigt = stress.voigt();

        // Final force-only evaluations need neither enhanced nor nodal Hessians
        if (!assemble_tangent) {
            const auto Bu       = this->green_lagrange_strain_displacement(enhanced_shape_derivatives, F);
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
        const Matrix6x13 Ba  = enhanced_green_lagrange_matrix(F, alpha_variations);
        const Matrix6x13 CBa = measure * C * Ba;

        // Assemble the enhanced residual and its material tangent. These blocks
        // are required by every local alpha Newton iteration.
        system.ra.noalias()  += measure * Ba.transpose() * stress_voigt;
        system.kaa.noalias() += Ba.transpose() * CBa;

        // The local Newton solve does not need any nodal tangent or residual block
        if (assemble_global_blocks) {
            const auto Bu                   = this->green_lagrange_strain_displacement(enhanced_shape_derivatives, F);
            const StaticMatrix<6, ndof> CBu = measure * C * Bu;
            system.ru.noalias()             += measure * Bu.transpose() * stress_voigt;
            system.kuu.noalias()            += Bu.transpose() * CBu;
            system.kua.noalias()            += Bu.transpose() * CBa;
            system.kau.noalias()            += Ba.transpose() * CBu;
        }

        // Material-only assembly is used when separating the condensed
        // geometric contribution from the complete finite-deformation tangent.
        if (!include_geometric) {
            continue;
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
                const Vec3 enhanced_dNb     = enhanced_shape_derivatives.row(node_b).transpose();
                const Precision coefficient = enhanced_dNa.dot(S * enhanced_dNb) * measure;

                for (Dim component = 0; component < D; ++component) {
                    system.kuu(D * node_a + component, D * node_b + component) += coefficient;
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
C3D8I::Vector13 C3D8I::solve_nonlinear_modes(
    const NonlinearPoints& points,
    const StaticVector<N>* thermal_free_strain
) {
    // Initialize the element-local enhanced state for the current trial geometry
    Vector13 alpha = Vector13::Zero();

    constexpr Index     max_iterations = 20;
    constexpr Precision tolerance      = Precision(1e-10);

    // Enforce stationarity of the thirteen local enhanced equations
    for (Index iteration = 0; iteration < max_iterations; ++iteration) {
        const EnhancedSystem system = assemble_nonlinear_system(
            points, alpha, false, false, true, true, thermal_free_strain);

        Eigen::FullPivLU<Matrix13> solver(system.kaa);
        logging::error(solver.isInvertible(),
            "C3D8I: singular local EAS tangent in element ", elem_id);

        // Apply the local Newton correction without exposing alpha globally
        const Vector13 delta = solver.solve(-system.ra);
        alpha                += delta;

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

/**
 * Collects element translations in node-major XYZ order.
 *
 * Connectivity selects the global field rows; the returned vector supplies the
 * local enhanced equations without changing model positions or material state.
 *
 * @param displacement Global nodal displacement field.
 * @return Element translational displacement vector.
 */
C3D8I::Vector24 C3D8I::local_displacement(const Field& displacement) {
    const StaticMatrix<N, D> local = this->nodal_data<D>(displacement);
    Vector24 result                = Vector24::Zero();

    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = local(node, dof);
        }
    }

    return result;
}

/**
 * Adds element translational forces to the global nodal accumulator.
 *
 * The output must use the NODE domain and provide at least XYZ components.
 * Each local node contributes according to its global connectivity identifier.
 *
 * @param node_forces Global nodal force accumulator.
 * @param local_force Element force in node-major XYZ ordering.
 */
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
 * Evaluates the C3D8I response through the common mechanical interface.
 *
 * All requests solve the thirteen local stationarity equations at the base
 * state u0 defined by linearization; nullptr denotes u0 = 0. The complete
 * condensed tangent is evaluated at u0. Internal force at a different requested
 * state u is returned through the affine approximation
 *
 *     f(u) = f(u0) + K_T(u0) (u - u0).
 *
 * The separately requested geometric stiffness is generated only by the
 * linearized PK2 stress increment from u0 to u. Existing stress at u0 remains
 * part of the complete base tangent. The perturbation geometric blocks are
 * condensed by differentiating the complete EAS Schur complement at u0.
 */
MapMatrix C3D8I::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force,
    const Field* displacement,
    const Field* linearization,
    const Field* thermal_free_strain,
    bool         update_state
) {
    const bool with_tangent   = tangent_buffer           != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force           != nullptr;

    if (!with_tangent && !with_geometric && !with_force) {
        return MapMatrix(nullptr, 0, 0);
    }

    logging::error(!with_force || displacement != nullptr,
        "C3D8I: internal force evaluation requires displacement");
    logging::error(!update_state || (linearization != nullptr && displacement == linearization),
        "C3D8I: material state requires an exact evaluation at the linearization state");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "C3D8I: thermal free strain must be scalar ELEMENT_NODAL data");

    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    Vector24 u_linearization = Vector24::Zero();
    if (linearization) u_linearization = local_displacement(*linearization);

    const Vector24 u = displacement ? local_displacement(*displacement) : u_linearization;
    const Vector24 delta = u - u_linearization;

    StaticMatrix<N, D> local_u;
    for (Index node = 0; node < N; ++node) {
        local_u.row(node) = u_linearization.template segment<3>(D * node).transpose();
    }

    const bool exact_thermal_state =
        thermal_free_strain && linearization != nullptr && displacement == linearization;

    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) =
                (*thermal_free_strain)(
                    static_cast<Index>(this->elem_nodal_offset) + node, 0);
        }
    }

    const StaticMatrix<N, D> current_coords = reference_coords + local_u;
    const NonlinearPoints points = nonlinear_points(reference_coords, current_coords);
    const Vector13 alpha = solve_nonlinear_modes(
        points,
        exact_thermal_state ? &nodal_thermal_strain : nullptr);

    const bool affine_force          = with_force && displacement != linearization;
    const bool need_complete_tangent = with_tangent || with_geometric || affine_force || thermal_free_strain;

    const EnhancedSystem system = assemble_nonlinear_system(
        points,
        alpha,
        update_state,
        true,
        need_complete_tangent,
        true,
        exact_thermal_state ? &nodal_thermal_strain : nullptr
    );

    Matrix24 complete = Matrix24::Zero();
    if (need_complete_tangent) {
        Eigen::FullPivLU<Matrix13> solver(system.kaa);
        logging::error(solver.isInvertible(),
            "C3D8I: singular converged local EAS tangent in element ", elem_id);
        complete = system.kuu - system.kua * solver.solve(system.kau);
    }

    if (with_geometric) {
        stiffness_geom(
            geometric_tangent_buffer,
            points,
            alpha,
            system,
            delta,
            exact_thermal_state ? nullptr : thermal_free_strain
        );
    }

    if (with_tangent) {
        MapMatrix mapped(tangent_buffer, ndof, ndof);
        mapped = complete;
    }

    if (with_force) {
        Vector24 force = system.ru;

        if (affine_force) {
            force.noalias() += complete * delta;
        }

        // Affine/reference paths retain the thermal source formulation.
        // Exact nonlinear paths already evaluated the constitutive law with
        // mechanical strain and must not subtract the same thermal state twice.
        if (thermal_free_strain && !exact_thermal_state) {
            const auto thermal = thermal_force(
                points, alpha, nodal_thermal_strain);
            Eigen::FullPivLU<Matrix13> solver(system.kaa);
            logging::error(solver.isInvertible(),
                "C3D8I: singular thermal condensation tangent in element ", elem_id);

            force -= thermal.first - system.kua * solver.solve(thermal.second);
        }

        assemble_local_force(*internal_force, force);
    }

    if (with_tangent) {
        return MapMatrix(tangent_buffer, ndof, ndof);
    }
    if (with_geometric) {
        return MapMatrix(geometric_tangent_buffer, ndof, ndof);
    }
    return MapMatrix(nullptr, 0, 0);
}

/**
 * Recovers exact or affine strain and Cauchy stress from the stationary EAS state.
 *
 * Exact recovery solves the enhanced stationarity equations at the requested
 * displacement. Affine recovery solves them at the expansion point and obtains
 * the enhanced increment from Kaa delta_alpha = -Kau delta_u. Reference thermal
 * loading adds its enhanced force source to this equation.
 *
 * Both paths use F_bar = F_c (I + sum alpha_m H_m), Green-Lagrange strain and PK2
 * constitutive stress. Affine recovery differentiates the complete Cauchy
 * push-forward, including det(F). Material history remains unchanged. Pointwise
 * results are copied to integration-point output or extrapolated to element nodes.
 *
 * @param strain Optional strain output field.
 * @param stress Optional Cauchy stress output field.
 * @param displacement Requested global nodal displacement.
 * @param rst Integration-point or natural nodal output coordinates.
 * @param offset First output row belonging to this element.
 * @param thermal_free_strain Optional scalar reference thermal free strain.
 * @param linearization Affine expansion displacement; null selects zero.
 */
void C3D8I::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    const Field*     linearization,
    const Field*     thermal_free_strain
) {
    const bool exact_state = linearization == &displacement;

    // Validate recovery coordinates and reference thermal loading
    logging::error(strain != nullptr || stress != nullptr,
        "C3D8I: stress/strain recovery requires at least one output field");
    logging::error(rst.cols() >= 3,
        "C3D8I: recovery coordinates require at least three columns");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "C3D8I: thermal free strain must be scalar ELEMENT_NODAL data");
    const auto& scheme = this->integration_scheme_stiffness();
    const RowMatrix ip_rst = this->stress_strain_ip_rst();
    const bool output_at_ip = rst.rows() == ip_rst.rows() && rst.leftCols(3).isApprox(ip_rst);
    logging::error(output_at_ip || rst.rows() == static_cast<Eigen::Index>(N),
        "C3D8I: stress/strain output must use integration points or element nodes");

    // Solve the same stationary finite enhanced state used by mechanical evaluation
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    const StaticMatrix<N, D> local_u = this->nodal_data<D>(displacement);
    StaticMatrix<N, D> local_state = StaticMatrix<N, D>::Zero();
    if (exact_state) {
        local_state = local_u;
    } else if (linearization) {
        local_state = this->nodal_data<D>(*linearization);
    }
    const StaticMatrix<N, D> local_delta = local_u - local_state;
    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) =
                (*thermal_free_strain)(
                    static_cast<Index>(this->elem_nodal_offset) + node, 0);
        }
    }

    const auto points = nonlinear_points(
        reference_coords,
        StaticMatrix<N, D>(reference_coords + local_state));
    const Vector13 alpha = solve_nonlinear_modes(
        points,
        (exact_state && thermal_free_strain) ? &nodal_thermal_strain : nullptr);

    // Differentiate local stationarity to recover the condensed enhanced increment
    Vector13 delta_alpha = Vector13::Zero();
    if (!exact_state) {
        const auto system = assemble_nonlinear_system(points, alpha, false, true, true);
        Vector24 delta_u;
        for (Index node = 0; node < N; ++node) delta_u.segment<D>(D * node) = local_delta.row(node).transpose();
        Vector13 residual = system.kau * delta_u;
        if (thermal_free_strain) residual -= thermal_force(points, alpha, nodal_thermal_strain).second;
        Eigen::FullPivLU<Matrix13> solver(system.kaa);
        logging::error(solver.isInvertible(),
            "C3D8I: singular enhanced recovery tangent in element ", elem_id);
        delta_alpha = -solver.solve(residual);
    }
    RowMatrix ip_strain = RowMatrix::Zero(scheme.count(), 6);
    RowMatrix ip_stress = RowMatrix::Zero(scheme.count(), 6);

    // Evaluate material response only at the exact state or affine expansion point
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto& point = points[ip];
        Mat3 enhancement = Mat3::Identity();
        Mat3 delta_enhancement = Mat3::Zero();
        for (Index mode = 0; mode < n_modes; ++mode) {
            enhancement       += alpha(mode) * point.modes[mode];
            delta_enhancement += delta_alpha(mode) * point.modes[mode];
        }
        const Mat3 F = point.compatible * enhancement;
        const auto green =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        Precision free = Precision(0);
        if (thermal_free_strain) {
            free = this->shape_function(
                point.natural(0), point.natural(1), point.natural(2))
                .dot(nodal_thermal_strain);
        }

        Vec6 constitutive_strain_values = green.voigt();
        if (exact_state && thermal_free_strain) {
            constitutive_strain_values.head<3>().array() -= free;
        }
        const VolumeStrainGreenLagrange constitutive_strain(
            constitutive_strain_values);

        const Precision* old_state =
            &(*this->_model_data->material_state_old)(this->mp_index(ip), 0);
        VolumeStressPK2 second_pk;
        Mat6 tangent;
        evaluate_material(
            point.natural(0), point.natural(1), point.natural(2),
            constitutive_strain, old_state, nullptr,
            second_pk, (exact_state && !thermal_free_strain) ? nullptr : &tangent);

        Vec6 thermal_stress = Vec6::Zero();
        if (thermal_free_strain && !exact_state) {
            thermal_stress =
                tangent * Vec6(free, free, free, 0, 0, 0);
        }

        const VolumeStressPK2 effective_pk(
            Vec6(second_pk.voigt() - thermal_stress));

        Vec6 recovered_strain = green.voigt();
        const Mat3 sigma = effective_pk.to_cauchy(F).tensor();
        Mat3 recovered_stress = sigma;
        if (!exact_state) {
            // Include compatible and stationary enhanced variations in delta F
            const Mat3 delta_compatible = local_delta.transpose() * point.derivatives;
            const Mat3 delta_F = delta_compatible * enhancement + point.compatible * delta_enhancement;
            const Mat3 delta_E = Precision(0.5) * (F.transpose() * delta_F + delta_F.transpose() * F);
            const Vec6 delta_strain = VolumeStrainGreenLagrange(delta_E).voigt();
            recovered_strain += delta_strain;
            const Vec6 mechanical_increment = delta_strain;
            const Mat3 S = effective_pk.tensor();
            const Mat3 delta_S =
                VolumeStressPK2(Vec6(tangent * mechanical_increment)).tensor();
            recovered_stress += (delta_F * S * F.transpose() + F * delta_S * F.transpose()
                               + F * S * delta_F.transpose()) / F.determinant()
                              - (F.inverse() * delta_F).trace() * sigma;
        }
        ip_strain.row(ip) = recovered_strain.transpose();
        ip_stress.row(ip) = VolumeStressCauchy(recovered_stress).voigt().transpose();
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
    const RowMatrix& E           = this->extrapolation_matrix();
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
