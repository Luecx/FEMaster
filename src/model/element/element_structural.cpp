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

} // namespace fem::model
