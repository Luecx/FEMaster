/**
 * @file truss.cpp
 * @brief Implements T3 topology, material access and elementary geometry.
 *
 * The two-node truss has only one geometric reference direction. Mechanical
 * state evaluation is implemented separately in truss_evaluate.cpp; mass and
 * distributed-field integration live in truss_integration.cpp; result recovery
 * lives in truss_recovery.cpp. This file therefore contains only the persistent
 * element interface and geometry that is shared by those operations.
 *
 * @see T3
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "truss.h"

namespace fem {
namespace model {

/**
 * Constructs a two-node truss from its element id and global node ids.
 *
 * Section, material, model data and material-point storage are bound later by
 * the compiled model infrastructure.
 *
 * @param elem_id Global element identifier.
 * @param node_ids_in Global identifiers of the two truss nodes.
 */
T3::T3(ID elem_id, std::array<ID, N> node_ids_in)
    : StructuralElement(elem_id),
      node_ids(node_ids_in) {}

ElDofs T3::dofs() const {
    return ElDofs{true, true, true, false, false, false};
}

Dim T3::dimensions() const {
    return 3;
}

Dim T3::n_nodes() const {
    return N;
}

Dim T3::num_ip() const {
    return 1;
}

const ID* T3::nodes() const {
    return node_ids.data();
}

SurfacePtr T3::surface(ID surface_id) {
    (void) surface_id;
    return nullptr;
}

std::string T3::type_name() const {
    return "T3";
}

/**
 * Resolves and validates the truss section assigned to the element.
 *
 * The section supplies the reference cross-sectional area A0 used by the
 * Total-Lagrangian virtual work
 *
 *     delta W_int = integral_0^L0 A0 S delta E dX.
 *
 * @return Assigned truss section.
 */
TrussSection* T3::get_section() {
    logging::error(this->_section != nullptr,
        "T3: missing section for element ", this->elem_id);

    auto* section = this->_section->template as<TrussSection>();
    logging::error(section != nullptr,
        "T3: section is not a truss section for element ", this->elem_id);
    return section;
}

/**
 * Resolves the material referenced by the assigned truss section.
 *
 * @return Material assigned through the truss section.
 */
material::MaterialPtr T3::get_material() {
    TrussSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "T3: no material set for element ", this->elem_id);
    return section->material_;
}

/**
 * Resolves the constitutive law required by the axial truss formulation.
 *
 * The returned material law must provide either the infinitesimal axial pair
 * (epsilon, sigma) or the finite Total-Lagrangian pair (E, S), depending on the
 * requested evaluation state.
 *
 * @return Non-owning pointer to the assigned elasticity model.
 */
material::Elasticity* T3::get_elasticity() {
    auto mat = get_material();
    logging::error(mat->has_elasticity(),
        "T3: material has no elasticity for element ", this->elem_id);
    return mat->elasticity().get();
}

/**
 * Returns one nodal position X_i from the undeformed reference configuration.
 *
 * @param local_node Local node index zero or one.
 * @return Reference position X_i.
 */
Vec3 T3::node_position_reference(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "T3: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "T3: reference positions field not set in model data");

    return this->_model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[local_node]));
}

/**
 * Returns one nodal position x_i from the current model configuration.
 *
 * This accessor is used only by geometry-based utilities such as current volume
 * and distributed-field integration. Mechanical evaluation does not infer its
 * state from this field; it constructs x0 = X + u0 explicitly from the supplied
 * linearization displacement.
 *
 * @param local_node Local node index zero or one.
 * @return Current model position x_i.
 */
Vec3 T3::node_position_current(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "T3: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
        "T3: current positions field not set in model data");

    return this->_model_data->positions->row_vec3(static_cast<Index>(node_ids[local_node]));
}

/**
 * Returns the undeformed truss length
 *
 *     L0 = ||X2 - X1||.
 *
 * @return Reference length L0.
 */
Precision T3::length_reference() const {
    return (node_position_reference(1) - node_position_reference(0)).norm();
}

/**
 * Returns the current geometric length
 *
 *     l = ||x2 - x1||.
 *
 * This quantity belongs to geometry/integration utilities. Mechanical
 * linearization computes its own base-state length from X + u0.
 *
 * @return Current geometric length l.
 */
Precision T3::length_current() const {
    return (node_position_current(1) - node_position_current(0)).norm();
}

/**
 * Returns the unit vector along the reference truss axis,
 *
 *     N0 = (X2 - X1) / L0.
 *
 * @return Reference axial direction N0.
 */
Vec3 T3::direction_reference() const {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in element ", this->elem_id);

    return (node_position_reference(1) - node_position_reference(0)) / L0;
}

/**
 * Returns the current geometric volume represented by the line element.
 *
 * The truss stores a reference area A0 and uses the current centerline length l,
 * so generic spatial integration employs
 *
 *     V = A0 l.
 *
 * This is intentionally distinct from the reference measure A0 L0 used by the
 * Total-Lagrangian mechanical weak form.
 *
 * @return Current geometric volume A0*l.
 */
Precision T3::volume() {
    return get_section()->area_ * length_current();
}

} // namespace model
} // namespace fem
