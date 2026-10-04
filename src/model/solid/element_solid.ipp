/**
 * @file element_solid.ipp
 * @brief Implements common solid geometry, material queries and thermal/mass operators.
 *
 * The template routines collect reference/current nodal geometry, construct
 * linearized and Total-Lagrangian strain-displacement operators, transform
 * constitutive response through `SolidSection` and integrate thermal and mass
 * matrices. Mechanical force and tangent assembly belongs to `evaluate()` in
 * `element_solid_compute.ipp`.
 *
 * State-neutral constitutive queries read globally enumerated committed
 * material-point rows and pass no persistent target state. Physical nonlinear
 * trial updates are performed only by `evaluate()`.
 *
 * @see SolidElement
 * @see SolidSection
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"
#include "../../section/section_solid.h"

namespace fem::model {

/**
 * Resolves the solid section assigned to this element.
 *
 * The section supplies material data and orientation for constitutive queries.
 * Missing assignments or sections of another type produce an error.
 *
 * @return Non-owning pointer to the assigned solid section.
 */
template<Index N>
SolidSection* SolidElement<N>::get_section() {
    // Validate the assigned section before exposing its constitutive data
    logging::error(this->_section != nullptr,
        "Section not set for element ", this->elem_id);
    auto section = this->_section->template as<SolidSection>();
    logging::error(section != nullptr,
        "Section is not a solid section for element ", this->elem_id);
    return section;
}

/**
 * Builds the additional element material rotation from its three Euler angles.
 *
 * The element-level orientation is composed with the section orientation during
 * constitutive evaluation. An absent orientation field returns the identity;
 * an existing field must contain three components.
 *
 * @return Rotation axes used in the reference material basis.
 */
template<Index N>
Mat3 SolidElement<N>::additional_material_rotation() const {
    // An absent element orientation leaves the section material basis unchanged
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return Mat3::Identity();
    }

    auto angles_field                       = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name, "': material orientation requires 3 components");

    const Vec3 angles = angles_field->row_vec3(static_cast<Index>(this->elem_id));
    return cos::RectangularSystem::euler(angles(0), angles(1), angles(2)).get_axes(Vec3::Zero());
}

/**
 * Gathers reference nodal positions in element connectivity order.
 *
 * The model must provide the reference-position field. These coordinates define
 * the undeformed geometry for material queries and Total-Lagrangian operators.
 *
 * @return One global reference XYZ coordinate row per element node.
 */
template<Index N>
auto SolidElement<N>::node_coords_reference() -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
        "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "reference positions field not set in model data");

    // Gather global positions in the element interpolation order
    const auto& positions = *this->_model_data->positions_reference;
    StaticMatrix<N, D> coords {};

    for (Index i = 0; i < N; ++i) {
        coords.row(i) = positions.row_vec3(static_cast<Index>(this->node_ids[i])).transpose();
    }

    return coords;
}

/**
 * Gathers current nodal positions in element connectivity order.
 *
 * The model must provide the current-position field. Volume, mass and distributed
 * loads use this model configuration rather than a supplied trial displacement.
 *
 * @return One global current XYZ coordinate row per element node.
 */
template<Index N>
auto SolidElement<N>::node_coords_current() -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
        "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
        "current positions field not set in model data");

    // Gather global positions in the element interpolation order
    const auto& positions = *this->_model_data->positions;
    StaticMatrix<N, D> coords {};

    for (Index i = 0; i < N; ++i) {
        coords.row(i) = positions.row_vec3(static_cast<Index>(this->node_ids[i])).transpose();
    }

    return coords;
}

/**
 * Reads the optional scalar stiffness factor assigned to this element.
 *
 * Constitutive stress and tangent use the same factor. Missing model data or an
 * absent scaling field returns one; an existing field must be scalar.
 *
 * @return Element topology stiffness multiplier.
 */
template<Index N>
Precision SolidElement<N>::element_stiffness_scale() const {
    // Elements without a topology scaling field retain their original response
    if (!this->_model_data || !this->_model_data->element_stiffness_scale) {
        return Precision(1);
    }

    auto scale_field                       = this->_model_data->element_stiffness_scale;
    logging::error(scale_field->components == 1,
        "Field '", scale_field->name, "': element stiffness scale requires 1 component");
    return (*scale_field)(static_cast<Index>(this->elem_id));
}

template<Index N>
Vec3 SolidElement<N>::material_position_reference(Precision r, Precision s, Precision t) {
    return this->interpolate<D>(this->node_coords_reference(), r, s, t);
}

/**
 * Constructs the infinitesimal strain-displacement operator from global gradients.
 *
 * Each node contributes three translational columns in node-major XYZ order.
 * Rows follow [epsilon_xx, epsilon_yy, epsilon_zz, gamma_yz, gamma_xz, gamma_xy],
 * with engineering shear strains.
 *
 * @param shape_der_global Spatial shape-function gradients, one row per node.
 * @return Six-by-3N matrix mapping nodal translations to engineering strain.
 */
template<Index N>
auto SolidElement<N>::strain_displacement(const StaticMatrix<N, D>& shape_der_global) -> StaticMatrix<n_strain, D * N> {
    // Prepare the six engineering-strain rows in node-major translation order
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();

    // Insert each node gradient into the normal and engineering-shear rows
    for (Index j = 0; j < N; j++) {
        Dim r1 = j * 3;
        Dim r2 = r1 + 1;
        Dim r3 = r1 + 2;

        B(0, r1) = shape_der_global(j, 0);
        B(1, r2) = shape_der_global(j, 1);
        B(2, r3) = shape_der_global(j, 2);

        B(3, r2) = shape_der_global(j, 2);
        B(3, r3) = shape_der_global(j, 1);

        B(4, r1) = shape_der_global(j, 2);
        B(4, r3) = shape_der_global(j, 0);

        B(5, r1) = shape_der_global(j, 1);
        B(5, r2) = shape_der_global(j, 0);
    }

    return B;
}

/**
 * Transforms natural shape derivatives into the supplied global geometry.
 *
 * With J(i,j) = dx_j/dxi_i, global derivatives are obtained from J^-1 dN/dxi.
 * Reference coordinates normally define dN/dX; the thermal-load path supplies
 * its current geometry to preserve its existing integration convention.
 *
 * @param reference_coords Global nodal geometry used for the transformation.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @param det Signed Jacobian determinant returned for the volume measure.
 * @param check_det Require a positive determinant when true. Recovery may disable
 * this check and reject invalid points itself.
 * @return N-by-three global shape-function gradients.
 */
template<Index N>
auto SolidElement<N>::shape_derivatives_reference(
    const StaticMatrix<N, D>& reference_coords,
    Precision                 r,
    Precision                 s,
    Precision                 t,
    Precision&                det,
    bool                      check_det)
    -> StaticMatrix<N, D> {
    // Evaluate the natural derivatives and the selected configuration Jacobian
    const StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    const StaticMatrix<D, D> J0              = jacobian(reference_coords, r, s, t);

    // Return the physical volume scaling and enforce positive orientation when requested
    det = J0.determinant();

    if (check_det) {
        logging::error(det > 0,
            "negative reference determinant encountered in element ", elem_id, "\ndet        : ", det,
            "\nCoordinates: ", reference_coords, "\nJacobi     : ", J0);
    }

    // Map the natural gradient covectors into global coordinates
    return (J0.inverse() * local_shape_der.transpose()).transpose();
}

/**
 * Builds the Total-Lagrangian derivative of Green-Lagrange strain.
 *
 * The reference gradients and deformation gradient form dE/du in node-major XYZ
 * ordering. Shear rows contain twice the tensor shear strain, matching the
 * work-conjugate engineering-Voigt material interface.
 *
 * @param dN_dX Shape-function gradients in global reference coordinates.
 * @param F Deformation gradient mapping reference vectors into current vectors.
 * @return Six-by-3N Green-Lagrange strain-displacement matrix.
 */
template<Index N>
auto SolidElement<N>::green_lagrange_strain_displacement(
    const StaticMatrix<N, D>& dN_dX,
    const Mat3&               F)
    -> StaticMatrix<n_strain, D * N> {
    // Prepare the six engineering-strain rows in node-major translation order
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();

    // Differentiate E = (F^T F - I)/2 with respect to each nodal translation
    for (Index a = 0; a < N; ++a) {
        for (Dim p = 0; p < D; ++p) {
            const Index col = D * a + p;

            B(0, col) = F(p, 0) * dN_dX(a, 0);
            B(1, col) = F(p, 1) * dN_dX(a, 1);
            B(2, col) = F(p, 2) * dN_dX(a, 2);

            B(3, col) = F(p, 1) * dN_dX(a, 2) + F(p, 2) * dN_dX(a, 1);
            B(4, col) = F(p, 2) * dN_dX(a, 0) + F(p, 0) * dN_dX(a, 2);
            B(5, col) = F(p, 0) * dN_dX(a, 1) + F(p, 1) * dN_dX(a, 0);
        }
    }

    return B;
}

/**
 * Evaluates the zero-strain Total-Lagrangian material tangent at one natural
 * coordinate.
 *
 * This helper is used for auxiliary constitutive stiffness queries such as
 * reduced-integration hourglass stabilization. Passing a null target state keeps
 * the evaluation state-neutral.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @return Global material tangent at zero Green-Lagrange strain.
 */
template<Index N>
auto SolidElement<N>::material_tangent_reference(
    Precision        r,
    Precision        s,
    Precision        t,
    const Precision* old_state,
    Precision*       new_state)
    -> StaticMatrix<n_strain, n_strain> {
    // Query the zero-strain PK2 tangent with the caller-selected history pointers
    VolumeStrainGreenLagrange zero_strain;
    VolumeStressPK2           zero_stress;
    Mat6                      tangent;
    evaluate_material(r, s, t, zero_strain, old_state, new_state, zero_stress, &tangent);
    return tangent;
}

/**
 * Interpolates K components of nodal data at a natural element coordinate.
 *
 * The topology-specific shape vector weights each component independently; data
 * rows must follow the same connectivity order as the shape functions.
 *
 * @tparam K Number of interpolated components.
 * @param data Element-local nodal values.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Interpolated component vector.
 */
template<Index N>
template<Dim K>
StaticVector<K> SolidElement<N>::interpolate(
    StaticMatrix<N, K> data,
    Precision          r,
    Precision          s,
    Precision          t) {
    // Weight each component with the topology interpolation at the natural point
    StaticMatrix<N, 1> shape_func = shape_function(r, s, t);
    StaticVector<K> res {};
    for (Index i = 0; i < K; i++) {
        res(i) = shape_func.dot(data.col(i));
    }
    return res;
}

/**
 * Builds the natural-to-global geometric Jacobian from nodal coordinates.
 *
 * Rows correspond to natural directions and columns to global XYZ directions:
 * J(i,j) = sum_a x_a,j dN_a/dxi_i. Its determinant supplies the signed physical
 * volume scaling and its inverse transforms natural shape derivatives.
 *
 * @param node_coords Global nodal coordinates in the selected configuration.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Three-by-three Jacobian with natural directions stored in rows.
 */
template<Index N>
auto SolidElement<N>::jacobian(
    const StaticMatrix<N, D>& node_coords,
    Precision                 r,
    Precision                 s,
    Precision                 t)
    -> StaticMatrix<D, D> {
    StaticMatrix<N, D> local_shape_derivative = shape_derivative(r, s, t);
    StaticMatrix<D, D> jacobian {};

    // Contract nodal coordinates with natural derivatives to form dx_j/dxi_i
    for (Dim m = 0; m < D; m++) {
        for (Dim n = 0; n < D; n++) {
            Precision dxn_drm = 0;
            for (Dim k = 0; k < N; k++) {
                dxn_drm += node_coords(k, n) * local_shape_derivative(k, m);
            }
            jacobian(m, n) = dxn_drm;
        }
    }

    return jacobian;
}

/**
 * Computes the deformation gradient from reference and current Jacobians.
 *
 * The Jacobian row convention gives F = J_current^T J_reference^-T. Positive
 * reference, current and deformation-gradient determinants are required before
 * the result is used for Green-Lagrange strain or stress push-forward.
 *
 * @param reference_coords Global nodal reference positions.
 * @param current_coords Global nodal positions at the evaluated configuration.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @return Global deformation gradient dx/dX.
 */
template<Index N>
Mat3 SolidElement<N>::deformation_gradient(
    const StaticMatrix<N, D>& reference_coords,
    const StaticMatrix<N, D>& current_coords,
    Precision                 r,
    Precision                 s,
    Precision                 t) {
    // Build both configuration mappings using the same natural point
    const Mat3 J_reference        = jacobian(reference_coords, r, s, t);
    const Mat3 J_current          = jacobian(current_coords, r, s, t);
    const Precision det_reference = J_reference.determinant();
    const Precision det_current   = J_current.determinant();

    logging::error(det_reference > Precision(0),
        "non-positive reference determinant encountered in element ", elem_id, "\ndet        : ", det_reference,
        "\nCoordinates: ", reference_coords, "\nJacobi     : ", J_reference);
    logging::error(det_current > Precision(0),
        "non-positive current determinant encountered in element ", elem_id, "\ndet        : ", det_current,
        "\nCoordinates: ", current_coords, "\nJacobi     : ", J_current);

    const Mat3 F          = J_current.transpose() * J_reference.inverse().transpose();
    const Precision det_F = F.determinant();
    logging::error(det_F > Precision(0),
        "non-positive deformation gradient determinant in element ", elem_id, "\ndet(F): ", det_F, "\nF     : ", F);
    return F;
}

/**
 * Collects selected global field components in element connectivity order.
 *
 * The field must contain the last requested component. Component j is read from
 * offset + j * stride; no field storage or material state is modified.
 *
 * @tparam K Number of components gathered per node.
 * @param full_data Global nodal field.
 * @param offset First selected component.
 * @param stride Component spacing.
 * @return N-by-K element-local nodal values.
 */
template<Index N>
template<Dim K>
StaticMatrix<N, K> SolidElement<N>::nodal_data(
    const Field& full_data,
    Index        offset,
    Index        stride) {
    StaticMatrix<N, K> res {};
    // Verify that every requested field component is available
    runtime_assert(
        full_data.components >= offset + stride * (K - 1) + 1,
        "cannot extract this many elements from the data"
    );

    // Gather global field rows in connectivity order with the selected component spacing
    for (Dim m = 0; m < N; m++) {
        for (Dim j = 0; j < K; j++) {
            Index n   = j * stride + offset;
            res(m, j) = full_data(static_cast<Index>(node_ids[m]), n);
        }
    }

    return res;
}

/**
 * Integrates the scalar-temperature conductivity operator in reference geometry.
 *
 * The material must provide conductivity k. Volume quadrature evaluates
 * K_T = integral (dN/dX) k (dN/dX)^T dV0, then removes numerical asymmetry and
 * copies the result into caller-owned storage. Constitutive history is unchanged.
 *
 * @param buffer Storage for the N-by-N thermal matrix.
 * @return Map onto the assembled conductivity matrix.
 */
template<Index N>
MapMatrix SolidElement<N>::conductivity(Precision* buffer) {
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    // Resolve the material conductivity used by the scalar-temperature operator
    auto* section = this->get_section();

    logging::error(section->material_->has_thermal_conductivity(),
        "Material has no thermal conductivity at element ", elem_id);

    auto cond = section->material_->get_thermal_conductivity();

    // Integrate conductivity against reference shape gradients over dV0
    std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &cond](Precision r, Precision s, Precision t) -> StaticMatrix<N, N> {
            Precision det0;
            const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(reference_coords, r, s, t, det0);
            return StaticMatrix<N, N>(dN_dX * cond * dN_dX.transpose() * det0);
    };
    StaticMatrix<N, N> conductivity = integration_scheme().integrate(func);

    // Remove numerical asymmetry from the analytically symmetric conductivity operator
    conductivity = 0.5 * (conductivity + conductivity.transpose());

    // Copy the integrated operator into caller-owned contiguous storage
    MapMatrix mapped{buffer, N, N};
    mapped = conductivity;
    return mapped;
}

/**
 * Integrates consistent thermal capacity in the reference configuration.
 *
 * Density rho and specific heat c_p must be assigned. Volume quadrature forms
 * C_T = integral rho c_p N N^T dV0. The symmetric scalar-temperature matrix
 * is written to caller-owned storage without modifying material history.
 *
 * @param buffer Storage for the N-by-N thermal matrix.
 * @return Map onto the assembled thermal capacity matrix.
 */
template<Index N>
MapMatrix SolidElement<N>::capacity(Precision* buffer) {
    // Collect reference geometry and validate density and heat capacity
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    auto* section = this->get_section();

    logging::error(section->material_->has_density(),
        "Material has no density at element ", elem_id);
    logging::error(section->material_->has_thermal_specific_heat(),
        "Material has no specific heat at element ", elem_id);

    const Precision rho = section->material_->get_density();
    const Precision cp  = section->material_->get_thermal_specific_heat();

    // Integrate rho c_p N N^T with the physical reference volume measure
    std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, rho, cp](Precision r, Precision s, Precision t) -> StaticMatrix<N, N> {
            Precision det0;

            const StaticMatrix<N, 1> Nf = this->shape_function(r, s, t);

            const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(reference_coords, r, s, t, det0);

            (void) dN_dX;

            return StaticMatrix<N, N>(Nf * Nf.transpose() * (rho * cp * det0));
    };

    StaticMatrix<N, N> capacity = integration_scheme().integrate(func);

    // Remove quadrature round-off asymmetry before copying to caller storage
    capacity = StaticMatrix<N, N>(0.5 * (capacity + capacity.transpose()));

    // Copy the integrated operator into caller-owned contiguous storage
    MapMatrix mapped{buffer, N, N};
    mapped = capacity;
    return mapped;
}

/**
 * Integrates the consistent translational mass matrix in current model geometry.
 *
 * The assigned material must provide density. Volume quadrature integrates the
 * scalar nodal products rho N_a N_b and expands them independently into the three
 * translational directions. Node-major XYZ ordering matches the mechanical DOFs.
 *
 * @param buffer Storage for the 3N-by-3N element matrix.
 * @return Map onto the consistent mass matrix.
 */
template<Index N>
MapMatrix SolidElement<N>::mass(Precision* buffer) {
    // Require material density before constructing the inertial operator
    logging::error(material() != nullptr,
        "no material assigned to element ", elem_id);
    logging::error(material()->has_density(),
        "material has no density assigned at element ", elem_id);

    Precision density = material()->get_density();

    // Integrate nodal shape products in the current model configuration
    StaticMatrix<N, D> node_coords = this->node_coords_current();

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, node_coords, density](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det;
            StaticMatrix<D, D> jac             = this->jacobian(node_coords, r, s, t);
            StaticMatrix<N, N> shape_func_mass = this->shape_function(r, s, t) * this->shape_function(r, s, t).transpose();
            det                                = jac.determinant();

            // Expand the mass matrix from N x N to D * N x D * N
            StaticMatrix<D * N, D * N> mass_local = StaticMatrix<D * N, D * N>::Zero();

            for (Index i = 0; i < N; i++) {
                for (Index j = 0; j < N; j++) {
                    for (Dim d = 0; d < D; d++) {
                        mass_local(D * i + d, D * j + d) = shape_func_mass(i, j);
                    }
                }
            }

            return mass_local * det * density;
    };

    StaticMatrix<D * N, D * N> mass = integration_scheme().integrate(func);

    // Copy the integrated operator into caller-owned contiguous storage
    MapMatrix mapped{buffer, D * N, D * N};
    mapped = mass;
    return mapped;
}

/**
 * Integrates signed physical volume in the current model configuration.
 *
 * The topology volume rule integrates det(J) over the natural element domain.
 * The signed measure is retained; this routine does not validate or repair an
 * inverted geometry.
 *
 * @return Current signed element volume.
 */
template<Index N>
Precision SolidElement<N>::volume() {
    // Use the current model geometry for the signed volume measure
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    std::function<Precision(Precision, Precision, Precision)> func =
        [this, node_coords_glob](Precision r, Precision s, Precision t) -> Precision {
            Precision det = jacobian(node_coords_glob, r, s, t).determinant();
            return det;
        };

    // Sum signed Jacobian determinants with the topology volume weights
    Precision volume = integration_scheme().integrate(func);
    return volume;
}
}  // namespace fem::model
