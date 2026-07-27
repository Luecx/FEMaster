//
// Created by f_eggers on 19.06.2026.
//
#include "element_structural.h"

#include <vector>

namespace fem::model {

MapMatrix StructuralElement::stiffness_tangent(Precision* buffer,
                                               Field&       ip_stress_state,
                                               NodeData&    nodal_forces,
                                               const Field& displacement) {
    const Index nip       = static_cast<Index>(num_ip());
    const int   ip_offset = static_cast<int>(ip_index(0));
    logging::error(nip > 0,
                   "Element ", elem_id,
                   ": tangent stiffness requires at least one integration point");

    compute_stress_state(ip_stress_state, displacement, ip_offset, true);

    MapMatrix tangent = stiffness(buffer);

    std::vector<Precision> geometric_storage(
        static_cast<size_t>(tangent.rows()) * static_cast<size_t>(tangent.cols()),
        Precision(0)
    );
    MapMatrix geometric = stiffness_geom(
        geometric_storage.data(),
        ip_stress_state,
        ip_offset
    );
    tangent += geometric;

    compute_internal_force_nonlinear(nodal_forces, ip_stress_state);

    return tangent;
}

/**
 * Evaluates and scatters the nonlinear internal force without assembling a
 * tangent matrix.
 *
 * The default implementation preserves the existing two-phase structural
 * interface: first the element stores its current nonlinear integration-point
 * stress or resultant state in the global `ip_stress_state` field, then the
 * internal force is recovered from that stored state. Elements with a combined
 * force-only material evaluation may override this function to avoid building
 * material and geometric tangent blocks during line-search trials.
 *
 * @param ip_stress_state Global integration-point stress/resultant field.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement defining the current configuration.
 */
void StructuralElement::internal_force_nonlinear(Field&       ip_stress_state,
                                                 NodeData&    nodal_forces,
                                                 const Field& displacement) {
    const Index nip       = static_cast<Index>(num_ip());
    const int   ip_offset = static_cast<int>(ip_index(0));

    logging::error(nip > 0,
        "Element ", elem_id,
        ": nonlinear internal force requires at least one integration point");

    compute_stress_state(ip_stress_state, displacement, ip_offset, true);
    compute_internal_force_nonlinear(nodal_forces, ip_stress_state);
}

} // namespace fem::model
