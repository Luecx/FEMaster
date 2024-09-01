/******************************************************************************
 * @file reduce_vec_to_vec.cpp
 * @brief Implements functions to reduce and expand vectors based on another vector.
 *
 * The `reduce` function removes elements from a vector `a` corresponding to
 * non-`NaN` elements in vector `b`, returning a reduced vector.
 *
 * The `expand` function expands a reduced vector back to the full size by
 * placing its values in the positions where vector `b` has `NaN` values, while
 * retaining non-`NaN` values from `b`.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#include "reduce_vec_to_vec.h"
#include "../core/core.h"
#include "assert.h"

namespace fem { namespace mattools {

DynamicVector reduce_vec_to_vec(const DynamicVector& a, const DynamicVector& b) {
    runtime_assert(a.size() == b.size(), "Input vectors must have the same size.");

    std::vector<double> reduced_values;

    // Iterate over both vectors
    for (size_t i = 0; i < a.size(); ++i) {
        // If the value in `b` is NaN, keep the corresponding value from `a`
        if (std::isnan(b(i))) {
            reduced_values.push_back(a(i));
        }
    }

    // Copy reduced values into a DynamicVector (Eigen::VectorXd)
    DynamicVector reduced_vector(reduced_values.size());
    for (int i = 0; i < reduced_values.size(); ++i) {
        reduced_vector(i) = reduced_values[i];
    }

    return reduced_vector;
}

DynamicVector expand_vec_to_vec(const DynamicVector& reduced_vector, const DynamicVector& b) {
    runtime_assert(reduced_vector.size() <= b.size(),"Reduced vector cannot be larger than vector b.");

    // Initialize the expanded vector to be the same size as vector `b`
    DynamicVector expanded_vector = b;

    int reduced_index = 0;

    // Iterate through the `b` vector
    for (size_t i = 0; i < b.size(); ++i) {
        // If the value in `b` is NaN, replace it with the next value from the reduced vector
        if (std::isnan(b(i)) && reduced_index < reduced_vector.size()) {
            expanded_vector(i) = reduced_vector(reduced_index);
            ++reduced_index;
        }
    }

    return expanded_vector;
}

} } // namespace fem::mattools
