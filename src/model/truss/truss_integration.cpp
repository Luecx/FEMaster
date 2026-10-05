/**
 * @file truss_integration.cpp
 * @brief Implements T3 mass, spatial integration and thermal load hooks.
 *
 * Mechanical stiffness is formulated on the reference measure A0 L0, whereas
 * generic spatial field integration uses the current geometric measure A0 l.
 * Keeping these operations separate makes that distinction explicit.
 *
 * @see T3
 *
 * @author Finn Eggers
 * @date 04.10.2026
 */

#include "truss.h"
#include <cmath>

namespace fem {
namespace model {
namespace {

/**
 * Returns the midpoint of the current truss axis,
 *
 *     x_m = 1/2 (x1 + x2).
 *
 * Distributed scalar, vector and tensor fields are sampled at this point. For
 * the current one-point truss integration this is the only spatial quadrature
 * location required.
 */
Vec3 midpoint(T3& element) {
    return Precision(0.5) *
        (element.node_position_current(0)
       + element.node_position_current(1));
}

/**
 * Returns the optional density multiplier used by generic field integration.
 *
 * Pure geometric integration uses one. Mass-weighted integration multiplies the
 * current volume measure by the assigned density rho.
 */
Precision density_scale(T3& element, bool scale_by_density) {
    if (!scale_by_density)
        return Precision(1);

    auto material = element.get_material();
    logging::error(material->has_density(),
       "T3: material density is required when scale_by_density=true for element ", element.elem_id);

    return material->get_density();
}

} // namespace

/**
 * Assembles the lumped translational mass matrix.
 *
 * The truss mass is evaluated from the reference volume
 *
 *     m = rho A0 L0.
 *
 * Half of the total mass belongs to each node. Because T3 has translational
 * degrees of freedom only, the same nodal mass is repeated on x, y and z:
 *
 *     M_i = (m/2) I_3.
 *
 * The reference measure is intentional: structural mass must not change merely
 * because the element stretches during a geometrically nonlinear analysis.
 *
 * @param buffer Caller-owned six-by-six matrix storage.
 * @return Map onto the lumped element mass matrix.
 */
MapMatrix T3::mass(Precision* buffer) {
    StaticMatrix<N * 3, N * 3> mass_matrix = StaticMatrix<N * 3, N * 3>::Zero();

    auto material = get_material();
    if (material->has_density()) {
        const Precision rho = material->get_density();
        const Precision A0  = get_section()->area_;
        const Precision L0  = length_reference();
        const Precision m   = rho * A0 * L0;

        for (Index node = 0; node < N; ++node) {
            mass_matrix.block(node * 3, node * 3, 3, 3) = Precision(0.5) * m * Mat3::Identity();
        }
    }

    MapMatrix mapped(buffer, N * 3, N * 3);
    mapped = mass_matrix;
    return mapped;
}

/**
 * Integrates a scalar field over the current truss volume.
 *
 * With one midpoint sample the spatial integral is approximated by
 *
 *     integral_V f(x) dV ~= f(x_m) A0 l.
 *
 * If density scaling is requested, the integrand is multiplied by rho:
 *
 *     integral_V rho f(x) dV ~= rho f(x_m) A0 l.
 */
Precision T3::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision A0 = get_section()->area_;
    const Precision l  = length_current();

    if (A0 <= Precision(0) || l <= Precision(0)) {
        return Precision(0);
    }

    return field(midpoint(*this))
         * density_scale(*this, scale_by_density)
         * A0
         * l;
}

/**
 * Integrates a vector field over the current truss volume,
 *
 *     integral_V f(x) dV ~= f(x_m) A0 l.
 *
 * Density scaling uses the same optional factor rho as the scalar integration
 * path.
 */
Vec3 T3::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision A0 = get_section()->area_;
    const Precision l  = length_current();

    if (A0 <= Precision(0) || l <= Precision(0)) {
        return Vec3::Zero();
    }

    return field(midpoint(*this))
         * density_scale(*this, scale_by_density)
         * A0
         * l;
}

/**
 * Integrates a distributed vector field and scatters its resultant equally to
 * the two nodes.
 *
 * The midpoint rule first forms
 *
 *     F = f(x_m) A0 l,
 *
 * optionally multiplied by rho. The two-node linear interpolation has equal
 * shape functions at the midpoint,
 *
 *     N1(0) = N2(0) = 1/2,
 *
 * hence each node receives F/2.
 */
void T3::integrate_vector_field(
    Field&          node_loads,
    bool            scale_by_density,
    const VecField& field
) {
    const Precision A0 = get_section()->area_;
    const Precision l  = length_current();

    if (A0 <= Precision(0) || l <= Precision(0)) {
        return;
    }

    const Vec3 force = field(midpoint(*this))
                     * density_scale(*this, scale_by_density)
                     * A0
                     * l;

    for (Index local_node = 0; local_node < N; ++local_node) {
        const Index node = static_cast<Index>(node_ids[local_node]);

        node_loads(node, 0) += Precision(0.5) * force(0);
        node_loads(node, 1) += Precision(0.5) * force(1);
        node_loads(node, 2) += Precision(0.5) * force(2);
    }
}

/**
 * Integrates a second-order tensor field over the current truss volume,
 *
 *     integral_V T(x) dV ~= T(x_m) A0 l.
 */
Mat3 T3::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision A0 = get_section()->area_;
    const Precision l  = length_current();

    if (A0 <= Precision(0) || l <= Precision(0)) {
        return Mat3::Zero();
    }

    return field(midpoint(*this))
         * density_scale(*this, scale_by_density)
         * A0
         * l;
}

/**
 * Assembles the positive equivalent nodal load generated by free axial
 * expansion of the current model temperature state.
 */
void T3::apply_thermal_expansion_load(Field& node_loads, const Field& node_temp) {
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
        "T3: thermal expansion requires a scalar NODE temperature field");
    logging::error(node_loads.domain == FieldDomain::NODE && node_loads.components >= 3,
        "T3: thermal expansion requires at least three nodal load components");

    auto material = get_material();
    if (!material->has_thermal_expansion()) {
        return;
    }

    const Precision zero = material->get_thermal_zero_temperature();
    const Precision alpha = material->get_thermal_expansion();

    Precision temperature = Precision(0);
    for (Index node = 0; node < N; ++node) {
        const Precision value = node_temp(static_cast<Index>(node_ids[node]), 0);
        temperature += std::isfinite(value) ? value : zero;
    }
    temperature /= static_cast<Precision>(N);

    const Precision free_strain = alpha * (temperature - zero);
    if (free_strain == Precision(0)) {
        return;
    }

    const Index state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    AxialStressPK2 stress;
    Precision tangent = Precision(0);
    get_elasticity()->evaluate(
        AxialStrainGreenLagrange(Precision(0)),
        old_state,
        nullptr,
        stress,
        &tangent
    );

    const Vec3 force = get_section()->area_ * tangent * free_strain * direction_reference();
    const Index node0 = static_cast<Index>(node_ids[0]);
    const Index node1 = static_cast<Index>(node_ids[1]);

    for (Dim d = 0; d < 3; ++d) {
        node_loads(node0, d) -= force(d);
        node_loads(node1, d) += force(d);
    }
}

/**
 * Stores the nodal scalar free thermal strain alpha (T - T0).
 */
void T3::apply_thermal_free_strain(Field& thermal_free_strain, const Field& node_temp) {
    logging::error(thermal_free_strain.domain == FieldDomain::ELEMENT_NODAL
                && thermal_free_strain.components == 1,
        "T3: thermal free strain requires scalar ELEMENT_NODAL storage");
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
        "T3: thermal free strain requires a scalar NODE temperature field");

    auto material = get_material();
    if (!material->has_thermal_expansion()) {
        return;
    }

    const Precision zero  = material->get_thermal_zero_temperature();
    const Precision alpha = material->get_thermal_expansion();

    for (Index node = 0; node < N; ++node) {
        const Precision value = node_temp(static_cast<Index>(node_ids[node]), 0);
        const Precision temperature = std::isfinite(value) ? value : zero;
        thermal_free_strain(static_cast<Index>(this->elem_nodal_offset) + node, 0) +=
            alpha * (temperature - zero);
    }
}

} // namespace model
} // namespace fem
