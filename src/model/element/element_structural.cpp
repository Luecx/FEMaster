/**
 * @file element_structural.cpp
 * @brief Implements non-inline functionality of the structural element base.
 *
 * `StructuralElement` defines the solver-facing mechanical interface shared by
 * structural formulations. The mechanical operators themselves depend on the
 * kinematics, integration rule and constitutive representation of each element
 * formulation and are therefore implemented by the derived element classes.
 *
 * In particular, nonlinear tangent and internal-force evaluation remain one
 * formulation-specific operation so that stress, constitutive tangent and other
 * integration-point quantities can stay local to the element evaluation.
 *
 * @see StructuralElement
 *
 * @author Finn Eggers
 * @date 01.09.2026
 */

#include "element_structural.h"

namespace fem::model {

// The structural base currently contains no out-of-line implementation. This
// translation unit remains the implementation counterpart of the interface and
// provides a natural location for functionality that is genuinely common to all
// structural formulations.

} // namespace fem::model
