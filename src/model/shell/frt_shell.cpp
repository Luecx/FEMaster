/**
 * @file frt_shell.cpp
 * @brief Builds and explicitly instantiates the common finite-rotation shell.
 *
 * The template implementation is divided into several implementation fragments
 * according to mathematical responsibility:
 *
 * - reference configuration and curved-surface geometry,
 * - finite-rotation kinematics and exact strain derivatives,
 * - shell-section response and consistent assembly,
 * - stress and strain result recovery,
 * - mass, volume and distributed-field integration.
 *
 * These files are included into this translation unit so that the complete
 * `FRTShell<N>` definition is visible during explicit instantiation. Only this
 * file needs to be added to the build for the common template implementation;
 * the included implementation fragments must not be compiled as independent
 * translation units.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell_reference.inl"
#include "frt_shell_kinematics.inl"
#include "frt_shell_assembly.inl"
#include "frt_shell_output.inl"
#include "frt_shell_integrate.inl"

namespace fem::model {

// Generate the common shell implementation for all supported topologies
template struct FRTShell<3>;
template struct FRTShell<4>;
template struct FRTShell<6>;
template struct FRTShell<8>;

} // namespace fem::model
