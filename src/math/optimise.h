#pragma once

#include "../core/types.h"
#include "quadrature.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace fem {
namespace optimise {

/**
 * @brief Perform optimization to find the closest point in isoparametric coordinates (r, s).
 *
 * This function finds the (r, s) coordinates on an isoparametric element surface that are closest
 * to a given 3D point. Depending on the domain type, the optimization is either constrained
 * or unconstrained. The optimization uses a simple gradient descent method and numerical
 * differentiation to find the optimal parameters.
 *
 * @param point_3d Target 3D point for which the closest isoparametric coordinates are to be found.
 * @param element_surface_func Lambda function that takes r, s and returns a 3D point on the surface.
 * @param domain ElementDomain specifying the type of the domain (DOMAIN_ISO_TRI or DOMAIN_ISO_QUAD).
 * @param constrained Boolean flag to specify if constraints should be applied (default: true).
 * @return Vector2d Optimized r, s coordinates in the isoparametric domain.
 */
StaticVector<2> find_isoparametric_point(StaticVector<3>                                      point_3d,
                                         std::function<StaticVector<3>(Precision, Precision)> surface_func,
                                         fem::quadrature::Domain                              domain,
                                         bool                                                 constrained = true);

}}  // namespace fem::optimise