/**
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
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
bool SolidElement<N>::has_material_orientation() const {
    const bool has_section_orientation =
        this->_section && this->_section->template as<SolidSection>() &&
        this->_section->template as<SolidSection>()->orientation_;
    const bool has_topo_orientation = this->_model_data && this->_model_data->material_orientation;
    return has_section_orientation || has_topo_orientation;
}

template<Index N>
Mat3 SolidElement<N>::section_orientation_basis(Precision r, Precision s, Precision t) {
    auto section = get_section();
    if (!section->orientation_) {
        return Mat3::Identity();
    }

    const Vec3 point_global = this->interpolate<D>(this->node_coords_reference(), r, s, t);
    const Vec3 point_local  = section->orientation_->to_local(point_global);
    return section->orientation_->get_axes(point_local);
}

template<Index N>
Mat3 SolidElement<N>::material_orientation_basis(Precision r, Precision s, Precision t) {
    Mat3 basis = section_orientation_basis(r, s, t);

    if (this->_model_data && this->_model_data->material_orientation) {
        auto angles_field = this->_model_data->material_orientation;
        logging::error(angles_field->components == 3,
                       "Field '", angles_field->name, "': material orientation requires 3 components");
        const Index row = static_cast<Index>(this->elem_id);
        const Vec3 angles = angles_field->row_vec3(row);
        const Mat3 topo_basis = cos::RectangularSystem::euler(angles(0), angles(1), angles(2)).get_axes(Vec3::Zero());
        basis = basis * topo_basis;
    }

    return basis;
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
Mat6 SolidElement<N>::material_strain_transformation(const Mat3& basis) {
    return VolumeStrain::get_transformation_matrix(Mat3::Identity(), basis);
}

template<Index N>
Mat6 SolidElement<N>::material_stress_transformation(const Mat3& basis) {
    return VolumeStress::get_transformation_matrix(basis, Mat3::Identity());
}

template<Index N>
Mat6 SolidElement<N>::material_strain_transformation_derivative(
    const Mat3& basis,
    const Mat3& basis_derivative
) {
    Mat6 derivative;

    for (Index component = 0; component < n_strain; ++component) {
        Vec6 unit = Vec6::Zero();
        unit(component) = Precision(1);

        const Mat3 global_tensor = VolumeStrain(unit).tensor();
        const Mat3 local_derivative =
            basis_derivative.transpose() * global_tensor * basis
            + basis.transpose() * global_tensor * basis_derivative;
        derivative.col(component) = VolumeStrain(local_derivative).voigt();
    }
    return derivative;
}

template<Index N>
Mat6 SolidElement<N>::material_stress_transformation_derivative(
    const Mat3& basis,
    const Mat3& basis_derivative
) {
    Mat6 derivative;

    for (Index component = 0; component < n_strain; ++component) {
        Vec6 unit = Vec6::Zero();
        unit(component) = Precision(1);

        const Mat3 local_tensor = VolumeStress(unit).tensor();
        const Mat3 global_derivative =
            basis_derivative * local_tensor * basis.transpose()
            + basis * local_tensor * basis_derivative.transpose();
        derivative.col(component) = VolumeStress(global_derivative).voigt();
    }
    return derivative;
}

template<Index N>
void SolidElement<N>::evaluate_material(Precision                     r,
                                        Precision                     s,
                                        Precision                     t,
                                        const VolumeStrainLinearized& global_strain,
                                        VolumeStressCauchy&           global_stress,
                                        Mat6&                         global_tangent) {
    logging::error(material() != nullptr, "no material assigned to element ", elem_id);
    logging::error(material()->has_elasticity(), "material has no elasticity assigned at element ", elem_id);

    auto elasticity = material()->elasticity();
    logging::error(elasticity->supports_volume_linearized(),
                   "material does not support linearized volume evaluation at element ", elem_id);
    logging::error(elasticity->state_size() == 0,
                   "SolidElement does not yet provide integration-point material state storage");

    const Mat3 basis = has_material_orientation()
        ? material_orientation_basis(r, s, t)
        : Mat3::Identity();

    const VolumeStrainLinearized local_strain =
        global_strain.transformed(Mat3::Identity(), basis);

    VolumeStressCauchy local_stress;
    Mat6               local_tangent;
    elasticity->evaluate(local_strain, nullptr, nullptr, local_stress, local_tangent);

    global_stress = local_stress.transformed(basis, Mat3::Identity());
    global_tangent = material_stress_transformation(basis)
        * local_tangent
        * material_strain_transformation(basis);

    if (this->_model_data && this->_model_data->element_stiffness_scale) {
        auto scale_field = this->_model_data->element_stiffness_scale;
        logging::error(scale_field->components == 1,
                       "Field '", scale_field->name, "': element stiffness scale requires 1 component");
        const Precision scaling = (*scale_field)(static_cast<Index>(this->elem_id));
        global_stress.voigt() *= scaling;
        global_tangent        *= scaling;
    }
}

template<Index N>
void SolidElement<N>::evaluate_material(Precision                        r,
                                        Precision                        s,
                                        Precision                        t,
                                        const VolumeStrainGreenLagrange& global_strain,
                                        VolumeStressPK2&                 global_stress,
                                        Mat6&                            global_tangent) {
    logging::error(material() != nullptr, "no material assigned to element ", elem_id);
    logging::error(material()->has_elasticity(), "material has no elasticity assigned at element ", elem_id);

    auto elasticity = material()->elasticity();
    logging::error(elasticity->supports_volume_green_lagrange(),
                   "material does not support Green-Lagrange volume evaluation at element ", elem_id);
    logging::error(elasticity->state_size() == 0,
                   "SolidElement does not yet provide integration-point material state storage");

    const Mat3 basis = has_material_orientation()
        ? material_orientation_basis(r, s, t)
        : Mat3::Identity();

    const VolumeStrainGreenLagrange local_strain =
        global_strain.transformed(Mat3::Identity(), basis);

    VolumeStressPK2 local_stress;
    Mat6            local_tangent;
    elasticity->evaluate(local_strain, nullptr, nullptr, local_stress, local_tangent);

    global_stress = local_stress.transformed(basis, Mat3::Identity());
    global_tangent = material_stress_transformation(basis)
        * local_tangent
        * material_strain_transformation(basis);

    if (this->_model_data && this->_model_data->element_stiffness_scale) {
        auto scale_field = this->_model_data->element_stiffness_scale;
        logging::error(scale_field->components == 1,
                       "Field '", scale_field->name, "': element stiffness scale requires 1 component");
        const Precision scaling = (*scale_field)(static_cast<Index>(this->elem_id));
        global_stress.voigt() *= scaling;
        global_tangent        *= scaling;
    }
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

template<Index N>
auto SolidElement<N>::material_tangent_reference(Precision r, Precision s, Precision t)
    -> StaticMatrix<n_strain, n_strain> {
    VolumeStrainGreenLagrange zero_strain;
    VolumeStressPK2           zero_stress;
    Mat6                      tangent;
    evaluate_material(r, s, t, zero_strain, zero_stress, tangent);
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

//-----------------------------------------------------------------------------
// stiffness
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::stiffness(Precision* buffer) {
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    StaticMatrix<N, D> current_coords = this->node_coords_current();

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &current_coords](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
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
            evaluate_material(r, s, t, strain, stress, C);
            StaticMatrix<D * N, D * N> res = B.transpose() * (C * B) * det0;
            return StaticMatrix<D * N, D * N>(res);
        };

    StaticMatrix<D * N, D * N> stiffness = integration_scheme().integrate(func);
    stiffness = 0.5 * (stiffness + stiffness.transpose()); // Symmetrize

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = stiffness;
    return mapped;
}

template<Index N>
MapMatrix
SolidElement<N>::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();

    Index ip_counter = 0;

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &reference_coords, &ip_stress, ip_start_idx, &ip_counter]
        (Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N>
    {
        const Index ip_row = static_cast<Index>(ip_start_idx) + ip_counter++;
        const VolumeStressPK2 stress{ip_stress.row_vec6(ip_row)};
        const Mat3 S = stress.tensor();

        Precision det0;
        const StaticMatrix<N, D> dN_dX =
            this->shape_derivatives_reference(reference_coords, r, s, t, det0);

        StaticMatrix<D * N, D * N> Kg = StaticMatrix<D * N, D * N>::Zero();

        for (Index a = 0; a < N; ++a) {
            Vec3 dNa;
            dNa << dN_dX(a, 0), dN_dX(a, 1), dN_dX(a, 2);

            for (Index b = 0; b < N; ++b) {
                Vec3 dNb;
                dNb << dN_dX(b, 0), dN_dX(b, 1), dN_dX(b, 2);

                const Precision s_ab = (dNa.transpose() * S * dNb)(0, 0);

                for (Dim d = 0; d < D; ++d) {
                    Kg(D * a + d, D * b + d) += s_ab * det0;
                }
            }
        }

        return Kg;
    };

    StaticMatrix<D * N, D * N> Kg = integration_scheme().integrate(func);
    Kg = 0.5 * (Kg + Kg.transpose()); // Symmetrize

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = Kg;
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
