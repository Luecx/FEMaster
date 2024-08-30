/******************************************************************************
 * @file assemble_bc.cpp
 * @brief Implementation of functions for handling and assembling constraints into NodeData,
 * with options for handling duplicate entries, either summing values or enforcing
 * consistency for constraints.
 *
 * This function implements the logic for handling duplicate entries based on the
 * provided options (ADD or SET). If inconsistencies are found, it throws an error.
 *
 * The function assumes that NodeData is defined as Eigen::Matrix<Precision, Eigen::Dynamic, 6>.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 ******************************************************************************/

#include "assemble_bc.h"

namespace fem { namespace mattools {

void assemble_bc(NodeData& target_bc, const NodeData& source_bc, DuplicateHandling handling) {
    if (target_bc.rows() != source_bc.rows() || target_bc.cols() != source_bc.cols()) {
        throw std::invalid_argument("BCMatrices must have the same dimensions for assembly.");
    }

    // Iterate through each entry of the NodeData
    for (int i = 0; i < target_bc.rows(); ++i) {
        for (int j = 0; j < target_bc.cols(); ++j) {
            const auto& source_val = source_bc(i, j);
            auto& target_val = target_bc(i, j);

            if (!std::isnan(source_val)) {  // Only consider non-NaN entries in source
                if (std::isnan(target_val)) {
                    // If target is NaN, simply assign the source value
                    target_val = source_val;
                } else {
                    // Handle based on the specified duplicate handling option
                    switch (handling) {
                        case DuplicateHandling::ADD:
                            // Add source value to the target
                            target_val += source_val;
                            break;
                        case DuplicateHandling::SET:
                            // Ensure consistency in constraints; use logging to check condition
                            logging::error(std::isnan(target_val) || target_val == source_val,
                                           "two different support collectors seem to both constrain the same "
                                           "node with different values, this is not allowed");
                            // If the values are inconsistent, logging will catch it, so we don't modify the value.
                            break;
                    }
                }
            }
        }
    }
}

} } // namespace fem::mattools
