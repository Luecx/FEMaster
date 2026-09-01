/**
 * @file reduce_mat_to_vec.cpp
 * @brief Implements conversion between nodal fields and active system vectors.
 *
 * The conversion is driven entirely by `SystemDofIds`. Negative entries mark
 * inactive nodal components, while every non-negative entry stores the exact
 * coefficient index in the compact algebraic vector. The implementation therefore
 * supports arbitrary node-by-component layouts, including the six-component
 * structural mapping and scalar node-by-one thermal systems.
 *
 * Reduction and expansion use the stored system identifiers rather than assuming
 * that active entries occur contiguously during field traversal. This preserves
 * the mapping semantics when system numbering is constructed or reordered by a
 * separate utility.
 *
 * @see numerate_dofs.h
 * @see SystemDofIds
 *
 * @author Finn Eggers
 * @date 28.08.2024
 */

#include "reduce_mat_to_vec.h"

#include "../core/logging.h"

namespace fem { namespace mattools {

/**
 * Reduces a nodal field into compact active-system ordering.
 *
 * The input field must cover the same nodal rows as the system mapping and expose
 * at least the mapped component count. The required vector size is obtained from
 * the largest non-negative system identifier. Every active mapping entry then
 * copies its corresponding physical field value to that exact vector location.
 *
 * Inactive entries are skipped and contribute no coefficient. Because the stored
 * identifiers are used directly, the function remains correct for arbitrary
 * valid system orderings and does not depend on row-major active-entry traversal.
 *
 * @param dof_ids Node-by-component global system identifiers with `-1` for
 *                inactive components.
 * @param field NODE-domain physical field containing the mapped components.
 * @return Compact vector indexed by the non-negative entries of `dof_ids`.
 */
DynamicVector reduce_mat_to_vec(const SystemDofIds& dof_ids, const model::Field& field) {
    // Validate that the physical field covers every node and component addressed
    // by the system mapping before any coefficient is copied.
    logging::error(field.domain == model::FieldDomain::NODE,
        "reduce_mat_to_vec: input field must use NODE domain");
    logging::error(field.rows == static_cast<Index>(dof_ids.rows()),
        "reduce_mat_to_vec: field row count does not match system DOF map");
    logging::error(field.components >= static_cast<Index>(dof_ids.cols()),
        "reduce_mat_to_vec: field has fewer components than system DOF map");

    // Allocate the compact algebraic vector from the highest active system index.
    // An empty mapping produces a zero-sized vector without a special traversal.
    const int system_size = dof_ids.size() == 0 ? 0 : dof_ids.maxCoeff() + 1;
    DynamicVector reduced = DynamicVector::Zero(system_size);

    // Traverse the complete nodal mapping and copy each active physical component
    // to the exact algebraic location stored in `dof_ids`.
    for (Eigen::Index row = 0; row < dof_ids.rows(); ++row) {
        for (Eigen::Index col = 0; col < dof_ids.cols(); ++col) {
            const int system_id = dof_ids(row, col);
            if (system_id >= 0) {
                reduced(system_id) = field(
                    static_cast<Index>(row),
                    static_cast<Index>(col)
                );
            }
        }
    }

    return reduced;
}

/**
 * Expands a compact active-system vector into nodal field storage.
 *
 * The returned NODE field has exactly the row and component dimensions of the
 * supplied system mapping and is initialized to zero. Every active mapping entry
 * copies the corresponding compact vector coefficient into its physical nodal
 * location, while inactive entries remain zero.
 *
 * The compact vector may contain additional coefficients, but it must cover the
 * complete active identifier range referenced by `dof_ids`.
 *
 * @param dof_ids Node-by-component global system identifiers with `-1` for
 *                inactive components.
 * @param reduced_vector Compact algebraic vector containing active coefficients.
 * @return NODE-domain field with dimensions matching `dof_ids`.
 */
model::Field expand_vec_to_mat(
    const SystemDofIds& dof_ids,
    const DynamicVector& reduced_vector
) {
    // Determine the largest coefficient referenced by the nodal mapping and
    // verify that the supplied compact vector covers the complete active range.
    const int required_size = dof_ids.size() == 0 ? 0 : dof_ids.maxCoeff() + 1;
    logging::error(reduced_vector.size() >= required_size,
        "expand_vec_to_mat: reduced vector is smaller than the active system DOF range");

    // Allocate nodal storage with exactly the mapping dimensions. Zero
    // initialization gives inactive components their defined expanded value.
    model::Field expanded{
        "EXPANDED_VECTOR",
        model::FieldDomain::NODE,
        static_cast<Index>(dof_ids.rows()),
        static_cast<Index>(dof_ids.cols())
    };
    expanded.set_zero();

    // Scatter every active compact coefficient to the nodal row and component
    // associated with its stored global system identifier.
    for (Eigen::Index row = 0; row < dof_ids.rows(); ++row) {
        for (Eigen::Index col = 0; col < dof_ids.cols(); ++col) {
            const int system_id = dof_ids(row, col);
            if (system_id >= 0) {
                expanded(
                    static_cast<Index>(row),
                    static_cast<Index>(col)
                ) = reduced_vector(system_id);
            }
        }
    }

    return expanded;
}

} } // namespace fem::mattools
