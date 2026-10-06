/**
 * @file element_solid_load.ipp
 * @brief Implements thermal loading and distributed-field integration for solids.
 *
 * Thermal expansion is converted into an equivalent nodal load by integrating
 * the reference material tangent against isotropic thermal strain. Auxiliary
 * load quadrature points reuse the closest constitutive integration-point state
 * row so every material call follows the same explicit state interface.
 *
 * The remaining routines integrate scalar, vector and tensor fields over the
 * current solid volume, optionally including material density.
 *
 * @see SolidElement
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

namespace fem::model {

/**
 * Integrates a scalar field over the current model volume.
 *
 * The topology volume rule evaluates the field at interpolated global positions
 * and multiplies by the signed current Jacobian determinant. Optional density
 * scaling requires an assigned material density. Constitutive history is not
 * evaluated or modified.
 *
 * @param scale_by_density Multiply the volume measure by material density.
 * @param field Callback evaluated at global current coordinates.
 * @return Volume integral of the supplied field.
 */
template<Index N>
Precision SolidElement<N>::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    // Gather current geometry for physical volume integration
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    // Apply material density only when requested by the caller
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = this->material();
        logging::error(mat != nullptr && mat->has_density(),
            "SolidElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    Precision result   = Precision(0);
    const auto& scheme = this->integration_scheme();
    // Evaluate field samples with the signed physical volume quadrature measure
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto pt     = scheme.get_point(ip);
        const Precision r = pt.r;
        const Precision s = pt.s;
        const Precision t = pt.t;
        const Precision w = pt.w;

        StaticMatrix<N, 1> Nvals   = this->shape_function(r, s, t);
        const StaticMatrix<D, D> J = this->jacobian(node_coords_glob, r, s, t);
        const Precision detJ       = J.determinant();

        // Interpolate the global sample position using current nodal coordinates
        Vec3 x_ip = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_ip += Nvals(i) * node_coords_glob.row(i);

        result += field(x_ip) * (rho * w * detJ);
    }
    return result;
}

/**
 * Integrates a vector field over the current model volume.
 *
 * The topology volume rule evaluates the field at interpolated global positions
 * and multiplies by the signed current Jacobian determinant. Optional density
 * scaling requires an assigned material density. Constitutive history is not
 * evaluated or modified.
 *
 * @param scale_by_density Multiply the volume measure by material density.
 * @param field Callback evaluated at global current coordinates.
 * @return Volume integral of the supplied field.
 */
template<Index N>
Vec3 SolidElement<N>::integrate_vector_field(bool scale_by_density, const VecField& field) {
    // Gather current geometry for physical volume integration
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    // Apply material density only when requested by the caller
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = this->material();
        logging::error(mat != nullptr && mat->has_density(),
            "SolidElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    Vec3 result        = Vec3::Zero();
    const auto& scheme = this->integration_scheme();
    // Evaluate field samples with the signed physical volume quadrature measure
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto pt     = scheme.get_point(ip);
        const Precision r = pt.r;
        const Precision s = pt.s;
        const Precision t = pt.t;
        const Precision w = pt.w;

        StaticMatrix<N, 1> Nvals   = this->shape_function(r, s, t);
        const StaticMatrix<D, D> J = this->jacobian(node_coords_glob, r, s, t);
        const Precision detJ       = J.determinant();

        // Interpolate the global sample position using current nodal coordinates
        Vec3 x_ip = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_ip += Nvals(i) * node_coords_glob.row(i);

        result += field(x_ip) * (rho * w * detJ);
    }
    return result;
}

/**
 * Integrates a distributed vector field into consistent nodal forces.
 *
 * The topology volume rule evaluates the field at interpolated global positions
 * and multiplies by the signed current Jacobian determinant. Optional density
 * scaling requires an assigned material density. Constitutive history is not
 * evaluated or modified.
 * Shape-function weights distribute each quadrature contribution to the global
 * nodal translational accumulator.
 *
 * @param node_loads Global nodal force accumulator.
 * @param scale_by_density Multiply the volume measure by material density.
 * @param field Callback evaluated at global current coordinates.
 */
template<Index N>
void SolidElement<N>::integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) {
    // Gather current geometry for physical volume integration
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    // Apply material density only when requested by the caller
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = this->material();
        logging::error(mat != nullptr && mat->has_density(),
            "SolidElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    const auto& scheme = this->integration_scheme();
    // Evaluate field samples with the signed physical volume quadrature measure
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto pt     = scheme.get_point(ip);
        const Precision r = pt.r;
        const Precision s = pt.s;
        const Precision t = pt.t;
        const Precision w = pt.w;

        StaticMatrix<N, 1> Nvals   = this->shape_function(r, s, t);
        const StaticMatrix<D, D> J = this->jacobian(node_coords_glob, r, s, t);
        const Precision detJ       = J.determinant();

        // Interpolate the global sample position using current nodal coordinates
        Vec3 x_ip = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_ip += Nvals(i) * node_coords_glob.row(i);

        Vec3 f_ip = field(x_ip) * (rho * w * detJ);

        for (Index i = 0; i < N; ++i) {
            const ID n_id       = this->node_ids[i];
            const Precision a   = Nvals(i);
            node_loads(n_id, 0) += a * f_ip(0);
            node_loads(n_id, 1) += a * f_ip(1);
            node_loads(n_id, 2) += a * f_ip(2);
        }
    }
}

/**
 * Integrates a tensor field over the current model volume.
 *
 * The topology volume rule evaluates the field at interpolated global positions
 * and multiplies by the signed current Jacobian determinant. Optional density
 * scaling requires an assigned material density. Constitutive history is not
 * evaluated or modified.
 *
 * @param scale_by_density Multiply the volume measure by material density.
 * @param field Callback evaluated at global current coordinates.
 * @return Volume integral of the supplied field.
 */
template<Index N>
Mat3 SolidElement<N>::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    // Gather current geometry for physical volume integration
    StaticMatrix<N, D> node_coords_glob = this->node_coords_current();

    // Apply material density only when requested by the caller
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = this->material();
        logging::error(mat != nullptr && mat->has_density(),
            "SolidElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    Mat3 result        = Mat3::Zero();
    const auto& scheme = this->integration_scheme();
    // Evaluate field samples with the signed physical volume quadrature measure
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto pt     = scheme.get_point(ip);
        const Precision r = pt.r;
        const Precision s = pt.s;
        const Precision t = pt.t;
        const Precision w = pt.w;

        StaticMatrix<N, 1> Nvals   = this->shape_function(r, s, t);
        const StaticMatrix<D, D> J = this->jacobian(node_coords_glob, r, s, t);
        const Precision detJ       = J.determinant();

        // Interpolate the global sample position using current nodal coordinates
        Vec3 x_ip = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_ip += Nvals(i) * node_coords_glob.row(i);

        result += field(x_ip) * (rho * w * detJ);
    }
    return result;
}

}  // namespace fem::model
