#pragma once

#include "../core/types.h"

#include <functional>
#include <ostream>

namespace fem {
namespace math {
namespace interpolate {

/**
 * \enum InterpolationFunction
 * \brief Enumerates the available interpolation methods.
 */
enum InterpolationFunction {
    CONSTANT,
    LINEAR,
    BILINEAR,
    QUADRATIC,
    BILINQUAD,
    CUBIC,
};

/**
 * \brief Performs interpolation and returns the result at a specified point.
 *
 * \tparam F The interpolation function.
 * \param[in] xyz Interpolation nodes.
 * \param[in] values The values at the interpolation nodes.
 * \param[in] center The point where interpolation result is needed.
 * \param[out] r2_values Optional output for R^2 values.
 * \return Matrix with interpolated values.
 *
 * \example
 * NodeData nodes = {{1,2,3}, {2,3,4}};
 * NodeData vals = {{5}, {6}};
 * Vec3 point = {1.5, 2.5, 3.5};
 * DynamicMatrix result = interpolate<LINEAR>(nodes, vals, point);
 */
template<InterpolationFunction F>
DynamicMatrix interpolate(const NodeData&        xyz,
                          const NodeData&        values,
                          const Vec3& center,
                          DynamicVector*         r2_values = nullptr);

/**
 * @brief Computes interpolated values using various interpolation methods.
 *
 * Given a set of nodes (xyz) and their corresponding values, this function interpolates to find
 * the value at a specific point (center). The function decides the interpolation method to use based on
 * the number of nodes, the accuracy factor, and the max allowed interpolation method.
 *
 * @param[in] xyz  Matrix of nodes with each row being a 3D coordinate.
 * @param[in] values Matrix of values corresponding to the nodes in xyz.
 * @param[in] center 3D coordinate where the interpolation result is desired.
 * @param[out] r2_values Vector that returns the R-squared values for the interpolation fit.
 * @param[in] accuracy_factor Factor that scales the number of nodes to determine the complexity of interpolation method
 * (lower values prefer simpler methods). Values above 1 should be avoided.
 * @param[in] max_accuracy The maximum interpolation method allowed.
 * @return Matrix of interpolated values.
 *
 * @example
 * NodeData xyz = ...; // Get some node data
 * NodeData values = ...; // Corresponding values for the nodes
 * Vec3 center = {x, y, z}; // Some point in space
 * DynamicVector r2;
 * auto result = fem::math::interpolate::interpolate(xyz, values, center, &r2, 1.0, InterpolationFunction::BILINEAR);
 * // Result now contains interpolated values at the center point.
 */
DynamicMatrix interpolate(const NodeData&        xyz,
                          const NodeData&        values,
                          const Vec3& center,
                          DynamicVector*         r2_values       = nullptr,
                          float                  accuracy_factor = 1,
                          InterpolationFunction  max_accuracy    = InterpolationFunction::CUBIC);

class Interpolator {
    private:
    InterpolationFunction method;
    float accuracy;

    public:
    // Constructor
    explicit Interpolator(InterpolationFunction method_ = InterpolationFunction::QUADRATIC,
                 float accuracy_ = 1);

    // Setter for interpolation method
    void set_function(InterpolationFunction method_);

    // Getter for interpolation method
    InterpolationFunction get_function() const;

    // Setter for accuracy
    void set_accuracy(float accuracy_);

    // Getter for accuracy
    float get_accuracy() const;

    // Perform interpolation
    DynamicMatrix operator()(const NodeData&        xyz,
                             const NodeData&        values,
                             const Vec3& center,
                             DynamicVector*         r2_values = nullptr);
};


}    // namespace interpolate
}    // namespace math
}    // namespace fem