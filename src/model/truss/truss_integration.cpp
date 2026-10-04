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
    return Precision(0.5)
         * (element.node_position_current(0) + element.node_position_current(1));
}

/**
 * Returns the optional density multiplier used by generic field integration.
 *
 * Pure geometric integration uses one. Mass-weighted integration multiplies the
 * current volume measure by the assigned density rho.
 */
Precision density_scale(T3& element, bool scale_by_density) {
    if (!scale_by_density) {
        return Precision(1);
    }

    auto material = element.get_material();
    logging::error(material && material->has_density(),
        "T3: material density is required when scale_by_density=true for element ",
        element.elem_id);

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
            mass_matrix.block(node * 3, node * 3, 3, 3) =
                Precision(0.5) * m * Mat3::Identity();
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
 * Applies equivalent nodal loading from a prescribed temperature field.
 *
 * Thermal expansion is not implemented for T3 yet. The function deliberately
 * remains a no-op while satisfying the common StructuralElement interface.
 */
void T3::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void) node_loads;
    (void) node_temp;
    (void) ref_temp;
}

} // namespace model
} // namespace fem
