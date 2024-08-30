/******************************************************************************
 * @file reduce_vec_to_vec.h
 * @brief Provides functions to reduce and expand vectors based on another vector.
 *
 * The `reduce` function removes elements from a vector `a` corresponding to
 * non-`NaN` elements in vector `b`, returning a reduced vector.
 *
 * The `expand` function expands a reduced vector back to the full size by
 * placing its values in the positions where vector `b` has `NaN` values, while
 * retaining non-`NaN` values from `b`.
 *
 * @param a The input vector from which elements will be reduced.
 * @param b The input vector that controls which elements in `a` will be discarded.
 * @param reduced_vector The reduced vector that was produced by `reduce`.
 * @return DynamicVector The reduced or expanded vector depending on the function.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types.h"
#include <Eigen/Dense>
#include <cmath>  // For std::isnan

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Reduces vector `a` by discarding elements that correspond to
 * non-`NaN` values in vector `b`.
 *
 * For each element in `b`, if the value is `NaN`, the corresponding element
 * from `a` is kept. Otherwise, it is discarded. The result is a reduced
 * vector of selected elements from `a`.
 *
 * @param a The input vector from which elements will be reduced.
 * @param b The input vector that indicates which elements in `a` to discard
 * (non-`NaN` values).
 * @return DynamicVector The reduced vector containing only the elements from `a`
 * corresponding to `NaN` values in `b`.
 ******************************************************************************/
DynamicVector reduce_vec_to_vec(const DynamicVector& a, const DynamicVector& b);

/******************************************************************************
 * @brief Expands a reduced vector back to its original size based on the `NaN`
 * values in vector `b`. The `NaN` values in `b` indicate where to insert the
 * values from the reduced vector.
 *
 * @param reduced_vector The reduced vector that was produced by `reduce`.
 * @param b The vector that controls where to place the values from the reduced vector.
 * @return DynamicVector The expanded vector, with values from the reduced vector
 * inserted at positions corresponding to `NaN` values in `b`, and the original
 * non-`NaN` values in `b` retained.
 ******************************************************************************/
DynamicVector expand_vec_to_vec(const DynamicVector& reduced_vector, const DynamicVector& b);

} } // namespace fem::mattools
