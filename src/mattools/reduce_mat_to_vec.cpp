/**
 * @file reduce_mat_to_vec.cpp
 * @brief Implements conversion between nodal fields and active system vectors.
 *
 * Inactive entries are identified by negative system indices. Active entries are
 * addressed by their stored global system identifier, so the conversion works
 * for arbitrary node-by-component mappings such as structural node-by-six and
 * scalar thermal node-by-one systems.
 *
 * @date Created on 28.08.2024
 */

#include "reduce_mat_to_vec.h"

#include "../core/logging.h"

namespace fem { namespace mattools {

DynamicVector reduce_mat_to_vec(const SystemDofIds& dof_ids, const model::Field& field) {
    logging::error(field.domain == model::FieldDomain::NODE,
        "reduce_mat_to_vec: input field must use NODE domain");
    logging::error(field.rows == static_cast<Index>(dof_ids.rows()),
        "reduce_mat_to_vec: field row count does not match system DOF map");
    logging::error(field.components >= static_cast<Index>(dof_ids.cols()),
        "reduce_mat_to_vec: field has fewer components than system DOF map");

    const int system_size = dof_ids.size() == 0 ? 0 : dof_ids.maxCoeff() + 1;
    DynamicVector reduced = DynamicVector::Zero(system_size);

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

model::Field expand_vec_to_mat(
    const SystemDofIds& dof_ids,
    const DynamicVector& reduced_vector
) {
    const int required_size = dof_ids.size() == 0 ? 0 : dof_ids.maxCoeff() + 1;
    logging::error(reduced_vector.size() >= required_size,
        "expand_vec_to_mat: reduced vector is smaller than the active system DOF range");

    model::Field expanded{
        "EXPANDED_VECTOR",
        model::FieldDomain::NODE,
        static_cast<Index>(dof_ids.rows()),
        static_cast<Index>(dof_ids.cols())
    };
    expanded.set_zero();

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
