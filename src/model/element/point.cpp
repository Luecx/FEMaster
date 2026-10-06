/**
 * @file point.cpp
 * @brief Implements the zero-dimensional structural point element.
 *
 * PointElement represents concentrated mass, rotary inertia and grounded spring
 * properties attached to one structural node. The element has no geometric
 * measure, integration points or constitutive material history.
 *
 * Its mechanical response is strictly linear:
 *
 *     f_int(u) = K u,
 *
 * with a constant diagonal spring tangent K. Since the element develops no
 * stress-dependent stiffness, its geometric tangent is identically zero.
 *
 * Concentrated translational mass and rotary inertia are represented by a
 * diagonal 6x6 mass matrix. Density-scaled generic field integration interprets
 * the section mass as the complete concentrated mass measure instead of a
 * density integrated over a finite geometric volume.
 *
 * A PointElement may temporarily exist without an assigned section while model
 * topology is being constructed or compiled. Topological queries therefore do
 * not require physical properties. Mechanical evaluation, mass evaluation and
 * density-scaled integration, however, require a valid PointMassSection.
 *
 * @see point.h
 * @see PointMassSection
 * @see StructuralElement
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "point.h"

namespace fem::model {

/**
 * Constructs a point element from its persistent topology.
 *
 * Only the element identifier and one-node connectivity belong to the persistent
 * definition. Section assignment, compiled offsets and ModelData binding are
 * supplied separately by the model construction and compilation process.
 *
 * @param elem_id Element identifier.
 * @param nodes One-entry array containing the connected node.
 */
PointElement::PointElement(ID elem_id, std::array<ID, N> nodes)
    : StructuralElement(elem_id), node_ids(nodes) {}

/**
 * Creates an independent copy of the persistent point-element definition.
 *
 * Runtime bindings are intentionally not copied. In particular, section
 * assignment, compiled storage offsets and ModelData binding are reconstructed
 * when the copied element is inserted into the compiled model.
 *
 * @return Independent point element with identical id and connectivity.
 */
ElementPtr PointElement::copy() const {
    return std::make_shared<PointElement>(elem_id, node_ids);
}

/**
 * Determines the structural degrees of freedom required by the point element.
 *
 * A PointElement is allowed to exist without an assigned section while semantic
 * topology is being assembled or while the compiled element copy has not yet
 * received its section assignment. Such an element has no known physical
 * properties and therefore activates no structural DOFs.
 *
 * Once a PointMassSection is assigned, translational mass activates all three
 * translational DOFs because concentrated inertia may act independently in every
 * spatial direction. Individual translational spring constants additionally
 * activate their corresponding translational DOFs.
 *
 * Rotary inertia and rotational spring constants analogously activate the three
 * rotational DOFs independently.
 *
 * @return Six-component structural DOF activation mask.
 */
ElDofs PointElement::dofs() const {
    // A section-less point element is valid during model construction. Until
    // physical properties are assigned it contributes no active DOFs.
    if (!_section) {
        return ElDofs{false, false, false, false, false, false};
    }

    // An assigned point section must use the dedicated point-mass formulation.
    const auto* section = _section->as<PointMassSection>();
    logging::error(section != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Translational mass activates all translational DOFs. Individual spring,
    // rotary-inertia and rotary-spring properties activate their corresponding
    // additional components.
    return ElDofs{
        section->mass_ != Precision(0) || section->spring_constants_(0) != Precision(0),
        section->mass_ != Precision(0) || section->spring_constants_(1) != Precision(0),
        section->mass_ != Precision(0) || section->spring_constants_(2) != Precision(0),
        section->rotary_inertia_(0) != Precision(0) || section->rotary_spring_constants_(0) != Precision(0),
        section->rotary_inertia_(1) != Precision(0) || section->rotary_spring_constants_(1) != Precision(0),
        section->rotary_inertia_(2) != Precision(0) || section->rotary_spring_constants_(2) != Precision(0)
    };
}

/**
 * Returns the element connectivity.
 *
 * PointElement contains exactly one global node.
 *
 * @return Pointer to the one-entry connectivity array.
 */
const ID* PointElement::nodes() const {
    return node_ids.data();
}

/**
 * Returns the stable textual point-element identifier.
 *
 * @return Element type name `POINT`.
 */
std::string PointElement::type_name() const {
    return "POINT";
}

/**
 * Returns the geometric volume of the point element.
 *
 * Concentrated properties do not represent a finite geometric region. Their
 * mass measure is handled explicitly by density-scaled field integration rather
 * than through a non-zero geometric volume.
 *
 * @return Always zero.
 */
Precision PointElement::volume() {
    return Precision(0);
}

/**
 * Evaluates the mechanical response of the point element.
 *
 * Mechanical evaluation requires an assigned PointMassSection. Section-less
 * PointElements are valid during model construction but have no physical
 * response that can be evaluated.
 *
 * The element represents independent translational and rotational springs
 * attached directly to ground. Its exact internal force is
 *
 *     f_int(u) = K u,
 *
 * and its tangent is the constant matrix
 *
 *     K_T = K.
 *
 * This is fully compatible with the StructuralElement base-state contract,
 * because for any base displacement u0,
 *
 *     f_int(u0) + K (u-u0)
 *       = K u0 + K (u-u0)
 *       = K u.
 *
 * The base displacement therefore does not have to be evaluated explicitly.
 * Temperature has no influence on the spring law.
 *
 * PointElement contains no stress-dependent mechanics, hence
 *
 *     K_G = 0.
 *
 * It also owns no constitutive material state, so update_state has no effect.
 *
 * Each output pointer may independently be null. Only requested quantities are
 * evaluated.
 *
 * @param tangent Optional caller-owned storage for the complete 6x6 tangent.
 * @param geometric_tangent Optional caller-owned storage for the geometric
 *        tangent, which is always zero.
 * @param internal_force Optional global nodal internal-force accumulator.
 * @param target_displacement Requested displacement state. Required when
 *        internal_force is requested.
 * @param target_temperature Ignored; point springs have no thermal response.
 * @param base_displacement Ignored because the spring law is exactly linear.
 * @param base_temperature Ignored because the spring law has no thermal response.
 * @param update_state Ignored because no constitutive history exists.
 * @return Mapping of the requested tangent buffer, geometric-tangent buffer or
 *         an empty matrix if no matrix output was requested.
 */
MapMatrix PointElement::evaluate(
    Precision*   tangent,
    Precision*   geometric_tangent,
    NodeData*    internal_force,
    const Field* target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature,
    bool         update_state
) {
    // Mechanical evaluation requires complete physical point properties.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Point springs are linear, temperature-independent and carry no persistent
    // constitutive history. Their response is therefore independent of the base
    // state and update_state.
    (void) target_temperature;
    (void) base_displacement;
    (void) base_temperature;
    (void) update_state;

    // Assemble the exact internal spring force
    //
    //     f_int = K u.
    //
    // Translational and rotational springs act independently on the six
    // structural DOFs.
    if (internal_force != nullptr) {
        logging::error(target_displacement != nullptr,
            "PointElement: internal force evaluation requires displacement");

        const Index node = static_cast<Index>(node_ids[0]);

        for (Index dof = 0; dof < 3; ++dof) {
            (*internal_force)(node, dof    ) += section->spring_constants_       (dof) * (*target_displacement)(node, dof);
            (*internal_force)(node, dof + 3) += section->rotary_spring_constants_(dof) * (*target_displacement)(node, dof + 3);
        }
    }

    // A linear grounded spring has no stress-dependent geometric stiffness.
    if (geometric_tangent != nullptr) {
        MapMatrix geometric(geometric_tangent, 6, 6);
        geometric.setZero();
    }

    // The complete tangent is the constant diagonal spring stiffness.
    if (tangent != nullptr) {
        MapMatrix result(tangent, 6, 6);
        result.setZero();

        // Translational spring stiffnesses.
        result(0, 0) = section->spring_constants_(0);
        result(1, 1) = section->spring_constants_(1);
        result(2, 2) = section->spring_constants_(2);

        // Rotational spring stiffnesses.
        result(3, 3) = section->rotary_spring_constants_(0);
        result(4, 4) = section->rotary_spring_constants_(1);
        result(5, 5) = section->rotary_spring_constants_(2);

        return result;
    }

    // If only the geometric tangent was requested, return its mapped storage.
    if (geometric_tangent != nullptr) {
        return MapMatrix(geometric_tangent, 6, 6);
    }

    // Internal-force-only evaluations require no local matrix representation.
    return MapMatrix(nullptr, 0, 0);
}

/**
 * Builds the concentrated 6x6 mass matrix.
 *
 * Mass evaluation requires an assigned PointMassSection. Translational mass is
 * isotropic and contributes the same scalar mass to all three translational
 * DOFs. Rotary inertia is stored component-wise for the three rotational DOFs.
 *
 * The current point-mass model contains no translation-rotation coupling and no
 * products of inertia, so the complete matrix is diagonal.
 *
 * @param buffer Caller-owned storage for the 6x6 element mass matrix.
 * @return Mapping of the assembled mass matrix.
 */
MapMatrix PointElement::mass(Precision* buffer) {
    // Mass evaluation requires complete physical point properties.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Initialize the complete local mass matrix.
    MapMatrix result(buffer, 6, 6);
    result.setZero();

    // Concentrated translational mass acts equally in all spatial directions.
    result(0, 0) = section->mass_;
    result(1, 1) = section->mass_;
    result(2, 2) = section->mass_;

    // Rotary inertia contributes independently to each rotational DOF.
    result(3, 3) = section->rotary_inertia_(0);
    result(4, 4) = section->rotary_inertia_(1);
    result(5, 5) = section->rotary_inertia_(2);

    return result;
}

/**
 * Integrates a scalar field over the point element.
 *
 * A zero-dimensional element has no geometric integration measure. Consequently
 * an unscaled scalar integral is always zero and does not require an assigned
 * section.
 *
 * When density scaling is requested, the concentrated section mass replaces the
 * distributed density-volume measure:
 *
 *     integral_V rho f(x) dV = m f(x_point).
 *
 * Density-scaled integration therefore requires a valid PointMassSection.
 *
 * @param scale_by_density Whether the concentrated mass measure is requested.
 * @param field Scalar field evaluated at the connected point.
 * @return Integrated scalar value.
 */
Precision PointElement::integrate_scalar_field(
    bool               scale_by_density,
    const ScalarField& field
) {
    // Without density scaling a point has no geometric measure.
    if (!scale_by_density) {
        return Precision(0);
    }

    // Density scaling requires the concentrated mass stored by the section.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Evaluate the field exactly at the point and multiply by the complete mass.
    return section->mass_ * field(node_position(0));
}

/**
 * Integrates a vector field over the point element.
 *
 * A zero-dimensional element has no unscaled geometric integration measure.
 * With density scaling, the usual distributed integral reduces to
 *
 *     integral_V rho b(x) dV = m b(x_point).
 *
 * Density-scaled integration therefore requires a valid PointMassSection.
 *
 * @param scale_by_density Whether the concentrated mass measure is requested.
 * @param field Vector field evaluated at the connected point.
 * @return Integrated three-component vector.
 */
Vec3 PointElement::integrate_vector_field(
    bool            scale_by_density,
    const VecField& field
) {
    // Without density scaling a point has no geometric measure.
    if (!scale_by_density) {
        return Vec3::Zero();
    }

    // Density scaling requires the concentrated mass stored by the section.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Evaluate the supplied vector field exactly at the connected point.
    return section->mass_ * field(node_position(0));
}

/**
 * Integrates a vector field and assembles it directly into a global nodal field.
 *
 * A zero-dimensional element has no unscaled geometric integration measure.
 * With density scaling,
 *
 *     F = m b(x_point).
 *
 * Because PointElement contains exactly one node, the complete integrated
 * resultant belongs to that node. No shape-function interpolation or nodal
 * distribution is necessary.
 *
 * Density-scaled integration requires a valid PointMassSection.
 *
 * @param node_loads Global nodal field receiving the integrated vector.
 * @param scale_by_density Whether the concentrated mass measure is requested.
 * @param field Vector field evaluated at the connected point.
 */
void PointElement::integrate_vector_field(
    Field&          node_loads,
    bool            scale_by_density,
    const VecField& field
) {
    // Without density scaling a point has no geometric measure.
    if (!scale_by_density) {
        return;
    }

    // Density scaling requires the concentrated mass stored by the section.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Evaluate the complete concentrated resultant at the connected point.
    const Index node  = static_cast<Index>(node_ids[0]);
    const Vec3  value = section->mass_ * field(node_position(0));

    // The one-node topology requires no further distribution.
    node_loads(node, 0) += value(0);
    node_loads(node, 1) += value(1);
    node_loads(node, 2) += value(2);
}

/**
 * Integrates a tensor field over the point element.
 *
 * A zero-dimensional element has no unscaled geometric integration measure.
 * With density scaling the concentrated mass replaces the distributed
 * density-volume measure:
 *
 *     integral_V rho A(x) dV = m A(x_point).
 *
 * Density-scaled integration therefore requires a valid PointMassSection.
 *
 * @param scale_by_density Whether the concentrated mass measure is requested.
 * @param field Tensor field evaluated at the connected point.
 * @return Integrated 3x3 tensor.
 */
Mat3 PointElement::integrate_tensor_field(
    bool            scale_by_density,
    const TenField& field
) {
    // Without density scaling a point has no geometric measure.
    if (!scale_by_density) {
        return Mat3::Zero();
    }

    // Density scaling requires the concentrated mass stored by the section.
    const PointMassSection* section = nullptr;

    logging::error(_section != nullptr,
        "PointElement: no section assigned to element ", elem_id);
    logging::error((section = _section->as<PointMassSection>()) != nullptr,
        "PointElement: section is not a PointMassSection for element ", elem_id);

    // Evaluate the tensor field exactly at the connected point.
    return section->mass_ * field(node_position(0));
}

} // namespace fem::model