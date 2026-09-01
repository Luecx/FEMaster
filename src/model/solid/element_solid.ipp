/**
 * @file element_solid.ipp
 * @brief Implements common solid geometry, interpolation and stiffness assembly.
 *
 * The template routines collect reference/current nodal geometry, construct
 * linearized and Total-Lagrangian strain-displacement operators, transform
 * constitutive response through `SolidSection` and integrate element matrices.
 *
 * State-neutral constitutive queries read globally enumerated committed
 * material-point rows and pass no persistent target state. Physical nonlinear
 * trial updates are performed only by `stiffness_tangent()`.
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
template<Index N>
SolidSection* SolidElement<N>::get_section() {
    logging::error(this->_section != nullptr, "Section not set for element ", this->elem_id);
    auto section = this->_section->template as<SolidSection>();
    logging::error(section != nullptr, "Section is not a solid section for element ", this->elem_id);
    return section;
}

template<Index N>
Mat3 SolidElement<N>::additional_material_rotation() const {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return Mat3::Identity();
    }

    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
                   "Field '", angles_field->name, "': material orientation requires 3 components");

    const Vec3 angles = angles_field->row_vec3(static_cast<Index>(this->elem_id));
    return cos::RectangularSystem::euler(
        angles(0),
        angles(1),
        angles(2)
    ).get_axes(Vec3::Zero());
}

template<Index N>
auto SolidElement<N>::node_coords_reference()
    -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
                   "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
                   "reference positions field not set in model data");

    const auto& positions = *this->_model_data->positions_reference;
    StaticMatrix<N, D> coords {};

    for (Index i = 0; i < N; ++i) {
        coords.row(i) =
            positions.row_vec3(static_cast<Index>(this->node_ids[i])).transpose();
    }

    return coords;
}

template<Index N>
auto SolidElement<N>::node_coords_current()
    -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
                   "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
                   "current positions field not set in model data");

    const auto& positions = *this->_model_data->positions;
    StaticMatrix<N, D> coords {};

    for (Index i = 0; i < N; ++i) {
        coords.row(i) =
            positions.row_vec3(static_cast<Index>(this->node_ids[i])).transpose();
    }

    return coords;
}

template<Index N>
Precision SolidElement<N>::element_stiffness_scale() const {
    if (!this->_model_data || !this->_model_data->element_stiffness_scale) {
        return Precision(1);
    }

    auto scale_field = this->_model_data->element_stiffness_scale;
    logging::error(scale_field->components == 1,
                   "Field '", scale_field->name, "': element stiffness scale requires 1 component");
    return (*scale_field)(static_cast<Index>(this->elem_id));
}

template<Index N>
Vec3 SolidElement<N>::material_position_reference(Precision r, Precision s, Precision t) {
    return this->interpolate<D>(this->node_coords_reference(), r, s, t);
}

/**
 * Evaluates linearized solid constitutive response at one natural point.
 *
 * The physical reference position, optional additional material rotation and
 * caller-selected state pointers are forwarded to `SolidSection`. Returned
 * global Cauchy stress and tangent are then multiplied by the optional scalar
 * topology stiffness field associated with this element.
 *
 * Passing `new_state == nullptr` requests a state-neutral constitutive query.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @param global_strain Linearized strain in global coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @param global_stress Cauchy stress returned in global coordinates.
 * @param global_tangent Consistent global material tangent.
 */
template<Index N>
void SolidElement<N>::evaluate_material(Precision                     r,
                                        Precision                     s,
                                        Precision                     t,
                                        const VolumeStrainLinearized& global_strain,
                                        const Precision*              old_state,
                                        Precision*                    new_state,
                                        VolumeStressCauchy&           global_stress,
                                        Mat6&                         global_tangent) {
    get_section()->evaluate(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        global_strain,
        old_state,
        new_state,
        global_stress,
        global_tangent
    );

    const Precision scaling = element_stiffness_scale();
    global_stress.voigt() *= scaling;
    global_tangent        *= scaling;
}

/**
 * Evaluates Total-Lagrangian solid constitutive response at one natural point.
 *
 * The selected old/new state pointers are passed directly to the section and
 * constitutive law. Green-Lagrange strain, returned PK2 stress and `dS/dE` are
 * all expressed in global reference coordinates. Optional topology scaling is
 * applied to both output quantities after constitutive evaluation.
 *
 * Passing `new_state == nullptr` requests a state-neutral constitutive query.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @param global_strain Green-Lagrange strain in global reference coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @param global_stress PK2 stress returned in global reference coordinates.
 * @param global_tangent Consistent global material tangent `dS/dE`.
 */
template<Index N>
void SolidElement<N>::evaluate_material(Precision                        r,
                                        Precision                        s,
                                        Precision                        t,
                                        const VolumeStrainGreenLagrange& global_strain,
                                        const Precision*                 old_state,
                                        Precision*                       new_state,
                                        VolumeStressPK2&                 global_stress,
                                        Mat6&                            global_tangent) {
    get_section()->evaluate(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        global_strain,
        old_state,
        new_state,
        global_stress,
        global_tangent
    );

    const Precision scaling = element_stiffness_scale();
    global_stress.voigt() *= scaling;
    global_tangent        *= scaling;
}

//-----------------------------------------------------------------------------
// strain_displacement
//-----------------------------------------------------------------------------
template<Index N>
auto SolidElement<N>::strain_displacement(const StaticMatrix<N, D>& shape_der_global)
    -> StaticMatrix<n_strain, D * N> {
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();

    for (Index j = 0; j < N; j++) {
        Dim r1   = j * 3;
        Dim r2   = r1 + 1;
        Dim r3   = r1 + 2;

        B(0, r1) = shape_der_global(j, 0);    // dx/dr
        B(1, r2) = shape_der_global(j, 1);    // dy/ds
        B(2, r3) = shape_der_global(j, 2);    // dz/dt

        B(3, r2) = shape_der_global(j, 2);    // dy/dt
        B(3, r3) = shape_der_global(j, 1);    // dz/ds

        B(4, r1) = shape_der_global(j, 2);    // dz/dr
        B(4, r3) = shape_der_global(j, 0);    // dx/dt

        B(5, r1) = shape_der_global(j, 1);    // dx/ds
        B(5, r2) = shape_der_global(j, 0);    // dy/dt
    }

    return B;
}

template<Index N>
auto SolidElement<N>::shape_derivatives_reference(const StaticMatrix<N, D>& reference_coords,
                                                  Precision r,
                                                  Precision s,
                                                  Precision t,
                                                  Precision& det,
                                                  bool check_det)
    -> StaticMatrix<N, D> {
    const StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    const StaticMatrix<D, D> J0 = jacobian(reference_coords, r, s, t);

    det = J0.determinant();

    if (check_det) {
        logging::error(det > 0, "negative reference determinant encountered in element ", elem_id,
                       "\ndet        : ", det,
                       "\nCoordinates: ", reference_coords,
                       "\nJacobi     : ", J0);
    }

    return (J0.inverse() * local_shape_der.transpose()).transpose();
}

template<Index N>
auto SolidElement<N>::green_lagrange_strain_displacement(const StaticMatrix<N, D>& dN_dX,
                                                         const Mat3& F)
    -> StaticMatrix<n_strain, D * N> {
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();

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
auto SolidElement<N>::material_tangent_reference(Precision r, Precision s, Precision t,
                                                 const Precision* old_state, Precision* new_state)
    -> StaticMatrix<n_strain, n_strain> {
    VolumeStrainGreenLagrange zero_strain;
    VolumeStressPK2           zero_stress;
    Mat6                      tangent;
    evaluate_material(r, s, t, zero_strain, old_state, new_state, zero_stress, tangent);
    return tangent;
}

template<Index N>
template<Dim K>
StaticVector<K> SolidElement<N>::interpolate(StaticMatrix<N, K> data, Precision r, Precision s, Precision t) {
    StaticMatrix<N, 1> shape_func = shape_function(r, s, t);
    StaticVector<K> res {};
    for (Index i = 0; i < K; i++) {
        res(i) = shape_func.dot(data.col(i));
    }
    return res;
}

//-----------------------------------------------------------------------------
// strain_displacements
//-----------------------------------------------------------------------------
template<Index N>
auto SolidElement<N>::strain_displacements(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t, Precision& det, bool check_det)
    -> StaticMatrix<n_strain, D * N> {
    StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    StaticMatrix<D, D> jac = jacobian(node_coords, r, s, t);

    det = jac.determinant();
    StaticMatrix<D, D> inv = jac.inverse();

    if (check_det) {
        logging::error(det > 0, "negative determinant encountered in element ", elem_id,
                       "\ndet        : ", det,
                       "\nCoordinates: ", node_coords,
                       "\nJacobi     : ", jac);
    }

    StaticMatrix<N, D> global_shape_der = (inv * local_shape_der.transpose()).transpose();
    return strain_displacement(global_shape_der);
}

//-----------------------------------------------------------------------------
// jacobian
//-----------------------------------------------------------------------------
template<Index N>
auto SolidElement<N>::jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t)
    -> StaticMatrix<D, D> {
    StaticMatrix<N, D> local_shape_derivative = shape_derivative(r, s, t);
    StaticMatrix<D, D> jacobian {};

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

template<Index N>
Mat3 SolidElement<N>::deformation_gradient(const StaticMatrix<N, D>& reference_coords,
                                           const StaticMatrix<N, D>& current_coords,
                                           Precision r,
                                           Precision s,
                                           Precision t) {
    const Mat3 J_reference = jacobian(reference_coords, r, s, t);
    const Mat3 J_current = jacobian(current_coords, r, s, t);
    const Precision det_reference = J_reference.determinant();
    const Precision det_current = J_current.determinant();

    logging::error(det_reference > Precision(0),
                   "non-positive reference determinant encountered in element ", elem_id,
                   "\ndet        : ", det_reference,
                   "\nCoordinates: ", reference_coords,
                   "\nJacobi     : ", J_reference);
    logging::error(det_current > Precision(0),
                   "non-positive current determinant encountered in element ", elem_id,
                   "\ndet        : ", det_current,
                   "\nCoordinates: ", current_coords,
                   "\nJacobi     : ", J_current);

    const Mat3 F = J_current.transpose() * J_reference.inverse().transpose();
    const Precision det_F = F.determinant();
    logging::error(det_F > Precision(0),
                   "non-positive deformation gradient determinant in element ", elem_id,
                   "\ndet(F): ", det_F,
                   "\nF     : ", F);
    return F;
}

//-----------------------------------------------------------------------------
// nodal_data
//-----------------------------------------------------------------------------
template<Index N>
template<Dim K>
StaticMatrix<N, K>
SolidElement<N>::nodal_data(const Field& full_data, Index offset, Index stride) {
    StaticMatrix<N, K> res {};
    runtime_assert(
        full_data.components >= offset + stride * (K - 1) + 1,
        "cannot extract this many elements from the data"
    );

    for (Dim m = 0; m < N; m++) {
        for (Dim j = 0; j < K; j++) {
            Index n = j * stride + offset;
            res(m, j) = full_data(static_cast<Index>(node_ids[m]), n);
        }
    }

    return res;
}

template<Index N>
MapMatrix
SolidElement<N>::conductivity(Precision* buffer) {
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    // store the section / material properties
    auto* section = this->get_section();

    logging::error(section->material_->has_thermal_conductivity(),
        "Material has no thermal conductivity at element ", elem_id);

    auto cond = section->material_->get_thermal_conductivity();

    // Evaluate one constitutive state row for every stiffness quadrature point
    std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &cond](Precision r, Precision s, Precision t) -> StaticMatrix<N, N> {
            Precision det0;
            const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(reference_coords, r, s, t, det0);
            return StaticMatrix<N, N>(dN_dX * cond * dN_dX.transpose() * det0);
    };
    StaticMatrix<N, N> conductivity = integration_scheme().integrate(func);

    // Remove only numerical asymmetry from the analytically symmetric material tangent
    conductivity = 0.5 * (conductivity + conductivity.transpose()); // Symmetrize

    MapMatrix mapped{buffer, N, N};
    mapped = conductivity;
    return mapped;
}

template<Index N>
MapMatrix
SolidElement<N>::capacity(Precision* buffer) {
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    auto* section = this->get_section();

    logging::error(section->material_->has_density(),
        "Material has no density at element ", elem_id);
    logging::error(section->material_->has_thermal_specific_heat(),
        "Material has no specific heat at element ", elem_id);

    const Precision rho = section->material_->get_density();
    const Precision cp  = section->material_->get_thermal_specific_heat();

    std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, rho, cp](Precision r, Precision s, Precision t) -> StaticMatrix<N, N> {
            Precision det0;

            const StaticMatrix<N, 1> Nf =
                this->shape_function(r, s, t);

            const StaticMatrix<N, D> dN_dX =
                this->shape_derivatives_reference(reference_coords, r, s, t, det0);

            (void) dN_dX;

            return StaticMatrix<N, N>(Nf * Nf.transpose() * (rho * cp * det0));
    };

    StaticMatrix<N, N> capacity =
        integration_scheme().integrate(func);

    capacity = StaticMatrix<N, N>(
        0.5 * (capacity + capacity.transpose())
    );

    MapMatrix mapped{buffer, N, N};
    mapped = capacity;
    return mapped;
}

/**
 * DEBUG A/B: restore the pre-refactor Total-Lagrangian solid stiffness path.
 *
 * This intentionally mirrors the previous implementation so the
 * shell_solid_tie regression can distinguish a change in solid stiffness from
 * the structural API refactor around it. It is not intended as the final design.
 */
template<Index N>
MapMatrix
SolidElement<N>::stiffness(Precision* buffer) {
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    StaticMatrix<N, D> current_coords   = this->node_coords_current();

    Index ip = 0;

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &current_coords, &ip](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det0;
            const StaticMatrix<N, D> dN_dX =
                this->shape_derivatives_reference(reference_coords, r, s, t, det0);
            const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);
            const StaticMatrix<n_strain, D * N> B =
                this->green_lagrange_strain_displacement(dN_dX, F);
            const VolumeStrainGreenLagrange strain =
                VolumeStrainGreenLagrange::from_deformation_gradient(F);
            VolumeStressPK2 stress;
            Mat6            C;

            const Index      state_row = this->mp_index(ip++);
            const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
            Precision*       new_state = &(*this->_model_data->material_state_new)(state_row, 0);
            evaluate_material(r, s, t, strain, old_state, new_state, stress, C);

            StaticMatrix<D * N, D * N> res = B.transpose() * (C * B) * det0;
            return StaticMatrix<D * N, D * N>(res);
        };

    StaticMatrix<D * N, D * N> stiffness = integration_scheme_stiffness().integrate(func);

    // Remove only numerical asymmetry from the analytically symmetric tangent.
    stiffness = Precision(0.5) * (stiffness + stiffness.transpose());

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = stiffness;
    return mapped;
}

/**
 * Integrates geometric stiffness for a supplied linearized prestress state.
 *
 * Prestress is reconstructed locally from the supplied nodal displacement field
 * rather than read from a global integration-point stress scratch field. At each
 * quadrature point
 *
 *     eps = B u_e,
 *     sigma = sigma(eps, state_committed),
 *
 * and the Cauchy stress contributes
 *
 *     K_geo,ab = integral grad(N_a)^T sigma grad(N_b) I_3 dV0.
 *
 * The constitutive evaluation is state-neutral and therefore cannot overwrite
 * the physical nonlinear trial material state.
 *
 * @param buffer Caller-provided dense element-matrix storage.
 * @param displacement Global nodal displacement field defining the prestress.
 * @return Mapped symmetric geometric stiffness matrix.
 */
template<Index N>
MapMatrix
SolidElement<N>::stiffness_geom(Precision* buffer, const Field& displacement) {
    const StaticMatrix<N, D> reference_coords  = this->node_coords_reference();
    const StaticMatrix<N, D> local_displacement = this->nodal_data<D>(displacement);
    const StaticMatrix<D, N> local_disp_mat(local_displacement.transpose());
    const auto local_displacement_vec =
        Eigen::Map<const StaticVector<D * N>>(local_disp_mat.data(), D * N);

    Index ip = 0;

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &local_displacement_vec, &ip]
        (Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N>
    {
        Precision det0;
        const StaticMatrix<N, D> dN_dX =
            this->shape_derivatives_reference(reference_coords, r, s, t, det0);
        const StaticMatrix<n_strain, D * N> B = this->strain_displacement(dN_dX);

        // Reconstruct the small-strain prestress directly from the supplied
        // displacement state. This auxiliary evaluation must remain state-neutral.
        const Vec6 strain_values = B * local_displacement_vec;
        const VolumeStrainLinearized strain(strain_values);
        const Index      state_row = this->mp_index(ip++);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        VolumeStressCauchy stress;
        Mat6               material_tangent;
        evaluate_material(r, s, t, strain, old_state, nullptr, stress, material_tangent);

        const Mat3 sigma = stress.tensor();
        StaticMatrix<D * N, D * N> geometric = StaticMatrix<D * N, D * N>::Zero();

        // Each scalar grad(N_a)^T sigma grad(N_b) coefficient multiplies the
        // translational identity block for the corresponding node pair.
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();

            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision stress_coefficient = dNa.dot(sigma * dNb) * det0;

                for (Dim d = 0; d < D; ++d) {
                    geometric(D * a + d, D * b + d) += stress_coefficient;
                }
            }
        }

        return geometric;
    };

    StaticMatrix<D * N, D * N> geometric = integration_scheme_stiffness().integrate(func);
    geometric = Precision(0.5) * (geometric + geometric.transpose());

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = geometric;
    return mapped;
}

//-----------------------------------------------------------------------------
// mass
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::mass(Precision* buffer) {
    logging::error(material() != nullptr, "no material assigned to element ", elem_id);
    logging::error(material()->has_density(), "material has no density assigned at element ", elem_id);

    Precision density = material()->get_density();

    StaticMatrix<N, D> node_coords = this->node_coords_current();

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, node_coords, density](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det;
            StaticMatrix<D, D> jac = this->jacobian(node_coords, r, s, t);
            StaticMatrix<N, N> shape_func_mass = this->shape_function(r, s, t) * this->shape_function(r, s, t).transpose();
            det = jac.determinant();

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

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = mass;
    return mapped;
}

//-----------------------------------------------------------------------------
// volume
//-----------------------------------------------------------------------------
template<Index N>
Precision
SolidElement<N>::volume() {
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    std::function<Precision(Precision, Precision, Precision)> func =
        [this, node_coords_glob](Precision r, Precision s, Precision t) -> Precision {
            Precision det = jacobian(node_coords_glob, r, s, t).determinant();
            return det;
        };

    Precision volume = integration_scheme().integrate(func);
    return volume;
}
}  // namespace fem::model