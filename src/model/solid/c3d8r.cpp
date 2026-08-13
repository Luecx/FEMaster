/**
 * @file c3d8r.cpp
 * @brief Implements C3D8R with finite-strain physical hourglass stabilization.
 *
 * The implementation is split into focused internal include files. The
 * Belytschko-Bindeman reference modal stiffness is embedded in an objective
 * Total-Lagrangian potential whose residual and tangent are differentiated
 * consistently.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

#include "c3d8r.h"

#include <cmath>
#include <vector>

namespace fem::model {

#include "c3d8r_base.ipp"
#include "c3d8r_hourglass_material.ipp"
#include "c3d8r_hourglass_reference.ipp"
#include "c3d8r_hourglass_nonlinear.ipp"
#include "c3d8r_assembly.ipp"

} // namespace fem::model
