//
// Created by f_eggers on 19.06.2026.
//
#include "element_structural.h"

#include <vector>

namespace fem::model {

MapMatrix StructuralElement::stiffness_tangent(Precision* buffer,
                                               NodeData& nodal_forces,
                                               const Field& displacement) {
    const Index nip = static_cast<Index>(num_ip());
    logging::error(nip > 0,
                   "Element ", elem_id,
                   ": tangent stiffness requires at least one integration point");

    Field stress_state{"ELEMENT_STRESS_STATE", FieldDomain::ELEMENT_IP, nip, 8};
    stress_state.set_zero();
    compute_stress_state(stress_state, displacement, 0, true);

    MapMatrix tangent = stiffness(buffer);

    std::vector<Precision> geometric_storage(
        static_cast<size_t>(tangent.rows()) * static_cast<size_t>(tangent.cols()),
        Precision(0)
    );
    MapMatrix geometric = stiffness_geom(geometric_storage.data(), stress_state, 0);
    tangent += geometric;

    compute_internal_force_nonlinear(nodal_forces, stress_state, 0);

    return tangent;
}

} // namespace fem::model
