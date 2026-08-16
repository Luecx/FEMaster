/**
 * @file element_solid.ipp
 * @brief Implements common solid geometry and constitutive dispatch.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"
#include "../../section/section_solid.h"

namespace fem::model {

template<Index N>
SolidSection* SolidElement<N>::get_section() {
    logging::error(this->_section != nullptr,
        "Section not set for element ", this->elem_id);
    auto section = this->_section->template as<SolidSection>();
    logging::error(section != nullptr,
        "Section is not a solid section for element ", this->elem_id);
    return section;
}

template<Index N>
Mat3 SolidElement<N>::additional_material_rotation() const {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return Mat3::Identity();
    }
    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name,
        "': material orientation requires 3 components");
    const Vec3 angles = angles_field->row_vec3(static_cast<Index>(this->elem_id));
    return cos::RectangularSystem::euler(
        angles(0), angles(1), angles(2)).get_axes(Vec3::Zero());
}

template<Index N>
auto SolidElement<N>::node_coords_reference() -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
        "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "reference positions field not set in model data");
    const auto& positions = *this->_model_data->positions_reference;
    StaticMatrix<N, D> coords {};
    for (Index i = 0; i < N; ++i) {
        coords.row(i) = positions.row_vec3(
            static_cast<Index>(this->node_ids[i])).transpose();
    }
    return coords;
}

template<Index N>
auto SolidElement<N>::node_coords_current() -> StaticMatrix<N, D> {
    logging::error(this->_model_data != nullptr,
        "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
        "current positions field not set in model data");
    const auto& positions = *this->_model_data->positions;
    StaticMatrix<N, D> coords {};
    for (Index i = 0; i < N; ++i) {
        coords.row(i) = positions.row_vec3(
            static_cast<Index>(this->node_ids[i])).transpose();
    }
    return coords;
}

template<Index N>
Precision SolidElement<N>::element_stiffness_scale() const {
    if (!this->_model_data || !this->_model_data->element_stiffness_scale) {
        return Precision(1);
    }
    auto field = this->_model_data->element_stiffness_scale;
    logging::error(field->components == 1,
        "Field '", field->name,
        "': element stiffness scale requires 1 component");
    return (*field)(static_cast<Index>(this->elem_id));
}

template<Index N>
Vec3 SolidElement<N>::material_position_reference(Precision r,
                                                  Precision s,
                                                  Precision t) {
    return this->interpolate<D>(this->node_coords_reference(), r, s, t);
}

template<Index N>
const Precision* SolidElement<N>::material_state_old(Index ip) const {
    logging::error(this->_section && this->_section->material_,
        "SolidElement requires material assignment");
    const Index size = this->_section->material_->constitutive_law().state_size();
    if (size == 0) return nullptr;
    logging::error(this->_model_data && this->_model_data->material_state_old,
        "Stateful constitutive update requires committed material state");
    logging::error(this->_model_data->material_state_old->components >= size,
        "Committed material-state field is too narrow");
    return &(*this->_model_data->material_state_old)(this->mp_index(ip), 0);
}

template<Index N>
Precision* SolidElement<N>::material_state_new(Index ip) const {
    logging::error(this->_section && this->_section->material_,
        "SolidElement requires material assignment");
    const Index size = this->_section->material_->constitutive_law().state_size();
    if (size == 0) return nullptr;
    logging::error(this->_model_data && this->_model_data->material_state_new,
        "Stateful constitutive update requires trial material state");
    logging::error(this->_model_data->material_state_new->components >= size,
        "Trial material-state field is too narrow");
    return &(*this->_model_data->material_state_new)(this->mp_index(ip), 0);
}

template<Index N>
void SolidElement<N>::update_material(Precision r,
                                      Precision s,
                                      Precision t,
                                      const VolumeStrainLinearized& strain,
                                      const Precision* state_old,
                                      Precision* state_new,
                                      VolumeStressCauchy& stress,
                                      Mat6& tangent) {
    get_section()->update(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        strain,
        state_old,
        state_new,
        stress,
        tangent
    );
    const Precision scaling = element_stiffness_scale();
    stress.voigt() *= scaling;
    tangent *= scaling;
}

template<Index N>
void SolidElement<N>::update_material(Precision r,
                                      Precision s,
                                      Precision t,
                                      const VolumeStrainGreenLagrange& strain,
                                      const Precision* state_old,
                                      Precision* state_new,
                                      VolumeStressPK2& stress,
                                      Mat6& tangent) {
    get_section()->update(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        strain,
        state_old,
        state_new,
        stress,
        tangent
    );
    const Precision scaling = element_stiffness_scale();
    stress.voigt() *= scaling;
    tangent *= scaling;
}

template<Index N>
void SolidElement<N>::recover_material(Precision r,
                                       Precision s,
                                       Precision t,
                                       const VolumeStrainLinearized& strain,
                                       const Precision* state,
                                       VolumeStressCauchy& stress) {
    get_section()->recover(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        strain,
        state,
        stress
    );
    stress.voigt() *= element_stiffness_scale();
}

template<Index N>
void SolidElement<N>::recover_material(Precision r,
                                       Precision s,
                                       Precision t,
                                       const VolumeStrainGreenLagrange& strain,
                                       const Precision* state,
                                       VolumeStressPK2& stress) {
    get_section()->recover(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        strain,
        state,
        stress
    );
    stress.voigt() *= element_stiffness_scale();
}

template<Index N>
auto SolidElement<N>::material_tangent_reference(Precision r,
                                                 Precision s,
                                                 Precision t)
    -> StaticMatrix<n_strain, n_strain> {
    Mat6 tangent = get_section()->elastic_tangent_reference(
        material_position_reference(r, s, t),
        additional_material_rotation()
    );
    tangent *= element_stiffness_scale();
    return tangent;
}

template<Index N>
auto SolidElement<N>::strain_displacement(
    const StaticMatrix<N, D>& shape_der_global
) -> StaticMatrix<n_strain, D * N> {
    StaticMatrix<n_strain, D * N> B =
        StaticMatrix<n_strain, D * N>::Zero();
    for (Index j = 0; j < N; ++j) {
        const Index x = j * 3;
        const Index y = x + 1;
        const Index z = x + 2;
        B(0, x) = shape_der_global(j, 0);
        B(1, y) = shape_der_global(j, 1);
        B(2, z) = shape_der_global(j, 2);
        B(3, y) = shape_der_global(j, 2);
        B(3, z) = shape_der_global(j, 1);
        B(4, x) = shape_der_global(j, 2);
        B(4, z) = shape_der_global(j, 0);
        B(5, x) = shape_der_global(j, 1);
        B(5, y) = shape_der_global(j, 0);
    }
    return B;
}

template<Index N>
auto SolidElement<N>::shape_derivatives_reference(
    const StaticMatrix<N, D>& reference_coords,
    Precision r,
    Precision s,
    Precision t,
    Precision& det,
    bool check_det
) -> StaticMatrix<N, D> {
    const StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    const StaticMatrix<D, D> J0 = jacobian(reference_coords, r, s, t);
    det = J0.determinant();
    if (check_det) {
        logging::error(det > 0,
            "negative reference determinant encountered in element ", elem_id,
            "\ndet        : ", det,
            "\nCoordinates: ", reference_coords,
            "\nJacobi     : ", J0);
    }
    return (J0.inverse() * local_shape_der.transpose()).transpose();
}

template<Index N>
auto SolidElement<N>::green_lagrange_strain_displacement(
    const StaticMatrix<N, D>& dN_dX,
    const Mat3& F
) -> StaticMatrix<n_strain, D * N> {
    StaticMatrix<n_strain, D * N> B =
        StaticMatrix<n_strain, D * N>::Zero();
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
template<Dim K>
StaticVector<K> SolidElement<N>::interpolate(StaticMatrix<N, K> data,
                                             Precision r,
                                             Precision s,
                                             Precision t) {
    const StaticMatrix<N, 1> shape_func = shape_function(r, s, t);
    StaticVector<K> result {};
    for (Index i = 0; i < K; ++i) result(i) = shape_func.dot(data.col(i));
    return result;
}

template<Index N>
auto SolidElement<N>::strain_displacements(const StaticMatrix<N, D>& node_coords,
                                           Precision r,
                                           Precision s,
                                           Precision t,
                                           Precision& det,
                                           bool check_det)
    -> StaticMatrix<n_strain, D * N> {
    const StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    const StaticMatrix<D, D> jac = jacobian(node_coords, r, s, t);
    det = jac.determinant();
    if (check_det) {
        logging::error(det > 0,
            "negative determinant encountered in element ", elem_id,
            "\ndet        : ", det,
            "\nCoordinates: ", node_coords,
            "\nJacobi     : ", jac);
    }
    const StaticMatrix<N, D> global_shape_der =
        (jac.inverse() * local_shape_der.transpose()).transpose();
    return strain_displacement(global_shape_der);
}

template<Index N>
auto SolidElement<N>::jacobian(const StaticMatrix<N, D>& node_coords,
                               Precision r,
                               Precision s,
                               Precision t) -> StaticMatrix<D, D> {
    const StaticMatrix<N, D> local_derivative = shape_derivative(r, s, t);
    StaticMatrix<D, D> result {};
    for (Dim m = 0; m < D; ++m) {
        for (Dim n = 0; n < D; ++n) {
            Precision value = 0;
            for (Index k = 0; k < N; ++k) {
                value += node_coords(k, n) * local_derivative(k, m);
            }
            result(m, n) = value;
        }
    }
    return result;
}

template<Index N>
Mat3 SolidElement<N>::deformation_gradient(
    const StaticMatrix<N, D>& reference_coords,
    const StaticMatrix<N, D>& current_coords,
    Precision r,
    Precision s,
    Precision t
) {
    const Mat3 J_reference = jacobian(reference_coords, r, s, t);
    const Mat3 J_current = jacobian(current_coords, r, s, t);
    const Precision det_reference = J_reference.determinant();
    const Precision det_current = J_current.determinant();
    logging::error(det_reference > Precision(0),
        "non-positive reference determinant encountered in element ", elem_id);
    logging::error(det_current > Precision(0),
        "non-positive current determinant encountered in element ", elem_id);
    const Mat3 F = J_current.transpose() * J_reference.inverse().transpose();
    logging::error(F.determinant() > Precision(0),
        "non-positive deformation gradient determinant in element ", elem_id);
    return F;
}

template<Index N>
template<Dim K>
StaticMatrix<N, K> SolidElement<N>::nodal_data(const Field& full_data,
                                               Index offset,
                                               Index stride) {
    StaticMatrix<N, K> result {};
    runtime_assert(
        full_data.components >= offset + stride * (K - 1) + 1,
        "cannot extract this many elements from the data"
    );
    for (Index m = 0; m < N; ++m) {
        for (Index j = 0; j < K; ++j) {
            result(m, j) = full_data(
                static_cast<Index>(node_ids[m]),
                j * stride + offset
            );
        }
    }
    return result;
}

template<Index N>
MapMatrix SolidElement<N>::stiffness(Precision* buffer) {
    const auto reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();
    StaticMatrix<D * N, D * N> K =
        StaticMatrix<D * N, D * N>::Zero();

    // Linear analyses use the initial elastic backbone, even if a plastic law is
    // attached. Material nonlinearity is integrated only by stiffness_tangent().
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        Precision det0;
        const auto dN_dX = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto B = strain_displacement(dN_dX);
        const Mat6 C = material_tangent_reference(point.r, point.s, point.t);
        K.noalias() += point.w * det0 * B.transpose() * C * B;
    }

    K = Precision(0.5) * (K + K.transpose());
    MapMatrix mapped{buffer, D * N, D * N};
    mapped = K;
    return mapped;
}

template<Index N>
MapMatrix SolidElement<N>::stiffness_geom(Precision* buffer,
                                          const Field& ip_stress,
                                          int ip_start_idx) {
    const auto reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();
    StaticMatrix<D * N, D * N> Kg =
        StaticMatrix<D * N, D * N>::Zero();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const Index row = static_cast<Index>(ip_start_idx) + ip;
        const Mat3 S = VolumeStressPK2(ip_stress.row_vec6(row)).tensor();
        Precision det0;
        const auto dN_dX = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Precision measure = point.w * det0;
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();
            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision value = dNa.dot(S * dNb) * measure;
                for (Dim d = 0; d < D; ++d) {
                    Kg(D * a + d, D * b + d) += value;
                }
            }
        }
    }

    Kg = Precision(0.5) * (Kg + Kg.transpose());
    MapMatrix mapped{buffer, D * N, D * N};
    mapped = Kg;
    return mapped;
}

template<Index N>
MapMatrix SolidElement<N>::mass(Precision* buffer) {
    logging::error(material() != nullptr,
        "no material assigned to element ", elem_id);
    logging::error(material()->has_density(),
        "material has no density assigned at element ", elem_id);
    const Precision density = material()->get_density();
    const StaticMatrix<N, D> node_coords = this->node_coords_current();

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, node_coords, density](Precision r, Precision s, Precision t) {
            const Precision det = this->jacobian(node_coords, r, s, t).determinant();
            const StaticMatrix<N, N> shape_mass =
                this->shape_function(r, s, t) * this->shape_function(r, s, t).transpose();
            StaticMatrix<D * N, D * N> local =
                StaticMatrix<D * N, D * N>::Zero();
            for (Index i = 0; i < N; ++i) {
                for (Index j = 0; j < N; ++j) {
                    for (Dim d = 0; d < D; ++d) {
                        local(D * i + d, D * j + d) = shape_mass(i, j);
                    }
                }
            }
            return local * det * density;
        };

    StaticMatrix<D * N, D * N> result = integration_scheme().integrate(func);
    MapMatrix mapped{buffer, D * N, D * N};
    mapped = result;
    return mapped;
}

template<Index N>
Precision SolidElement<N>::volume() {
    const StaticMatrix<N, D> coords = this->node_coords_current();
    std::function<Precision(Precision, Precision, Precision)> func =
        [this, coords](Precision r, Precision s, Precision t) {
            return jacobian(coords, r, s, t).determinant();
        };
    return integration_scheme().integrate(func);
}

} // namespace fem::model
