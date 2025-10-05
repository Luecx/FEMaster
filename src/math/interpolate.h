/******************************************************************************
 * @file interpolate.h
 * @brief Declares interpolation helpers for scattered nodal data.
 *
 * The API provides compile-time polynomial interpolation routines and a small
 * wrapper class that selects the interpolation order based on configuration.
 *
 * @see src/math/interpolate.cpp
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace math {
namespace interpolate {

/******************************************************************************
 * @enum InterpolationFunction
 * @brief Enumerates supported interpolation polynomials ordered by complexity.
 ******************************************************************************/
enum InterpolationFunction {
    CONSTANT,
    LINEAR,
    BILINEAR,
    QUADRATIC,
    BILINQUAD,
    CUBIC,
};

/******************************************************************************
 * @brief Interpolates values at `center` using the selected polynomial order.
 *
 * @tparam F Interpolation function to evaluate.
 * @param xyz   Positions of the interpolation nodes.
 * @param values Values associated with the nodes.
 * @param center Evaluation point in 3D space.
 * @param r2_values Optional output containing R² statistics per column.
 * @return Interpolated values stored as a single-row matrix.
 ******************************************************************************/
template<InterpolationFunction F>
DynamicMatrix interpolate(const NodeData& xyz,
                          const NodeData& values,
                          const Vec3& center,
                          DynamicVector* r2_values = nullptr);

/******************************************************************************
 * @brief Chooses an interpolation order based on available nodes and settings.
 *
 * @param xyz   Positions of the interpolation nodes.
 * @param values Values associated with the nodes.
 * @param center Evaluation point in 3D space.
 * @param r2_values Optional output containing R² statistics per column.
 * @param accuracy_factor Scaling applied to the node count when selecting the
 *                         interpolation order (values > 1 favour higher orders).
 * @param max_accuracy Hard upper limit on the interpolation order.
 * @return Interpolated values stored as a single-row matrix.
 ******************************************************************************/
DynamicMatrix interpolate(const NodeData& xyz,
                          const NodeData& values,
                          const Vec3& center,
                          DynamicVector* r2_values = nullptr,
                          float accuracy_factor = 1.0F,
                          InterpolationFunction max_accuracy = InterpolationFunction::CUBIC);

/******************************************************************************
 * @class Interpolator
 * @brief Thin wrapper that stores interpolation settings for repeated use.
 ******************************************************************************/
class Interpolator {
public:
    /// Constructs an interpolator with the supplied method and accuracy factor.
    explicit Interpolator(InterpolationFunction method = InterpolationFunction::QUADRATIC,
                          float accuracy = 1.0F);

    /// Sets the interpolation function to use.
    void set_function(InterpolationFunction method);

    /// Returns the currently configured interpolation function.
    InterpolationFunction get_function() const;

    /// Sets the accuracy factor.
    void set_accuracy(float accuracy);

    /// Returns the configured accuracy factor.
    float get_accuracy() const;

    /******************************************************************************
     * @brief Executes the interpolation using the stored configuration.
     ******************************************************************************/
    DynamicMatrix operator()(const NodeData& xyz,
                             const NodeData& values,
                             const Vec3& center,
                             DynamicVector* r2_values = nullptr);

private:
    InterpolationFunction m_method;
    float m_accuracy;
};

} // namespace interpolate
} // namespace math
} // namespace fem

