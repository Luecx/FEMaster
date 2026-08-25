/**
 * @file model_point_damping.cpp
 * @brief Assembles viscous ground damping contributed by point elements.
 *
 * Point-ground damping is deliberately local to `PointMassSection` rather than
 * part of the generic structural-element interface. The current implementation
 * supports only diagonal translational and rotational damping against ground, so
 * assembling those coefficients directly keeps the broader element hierarchy
 * free of a damping callback that most formulations do not need.
 *
 * @see PointElement
 * @see PointMassSection
 * @see Model::build_point_damping_matrix
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "model.h"

#include "element/point.h"

namespace fem::model {

/**
 * Assembles the global viscous damping matrix of all compiled point elements.
 *
 * Each `PointElement` contributes the six diagonal coefficients stored in its
 * `PointMassSection`. Inactive generalized DOFs have negative global indices and
 * are skipped. Non-point elements do not participate; distributed/Rayleigh
 * damping remains a separate load-case-level contribution.
 *
 * @param indices Node-by-six mapping from generalized nodal DOFs to global
 *                system indices.
 * @return Sparse global damping matrix in the same active-DOF space as K and M.
 */
SparseMatrix Model::build_point_damping_matrix(SystemDofIds& indices) {
    const Index global_size = static_cast<Index>(indices.maxCoeff() + 1);
    TripletList triplets;

    // One point element contributes at most six diagonal entries.
    triplets.reserve(_data->elements.size() * 6);

    for (const auto& element : _data->elements) {
        if (!element) continue;

        const auto* point = element->as<PointElement>();
        if (!point || !point->_section) continue;

        const auto* section = point->_section->as<PointMassSection>();
        logging::error(section != nullptr,
            "Model: point element ", point->elem_id, " has a non-point-mass section");

        const ID node = point->nodes()[0];
        for (Index dof = 0; dof < 3; ++dof) {
            const Index translation = indices(node, dof);
            const Index rotation    = indices(node, dof + 3);

            if (translation >= 0 && section->damping_constants_(dof) != Precision(0)) {
                triplets.emplace_back(translation, translation, section->damping_constants_(dof));
            }
            if (rotation >= 0 && section->rotary_damping_constants_(dof) != Precision(0)) {
                triplets.emplace_back(rotation, rotation, section->rotary_damping_constants_(dof));
            }
        }
    }

    SparseMatrix damping(global_size, global_size);
    damping.setFromTriplets(triplets.begin(), triplets.end());
    damping.makeCompressed();
    return damping;
}

} // namespace fem::model
