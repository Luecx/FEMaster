/**
 * @file plasticity.cpp
 * @brief Implements default unsupported plastic constitutive operations.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "plasticity.h"

#include "../core/logging.h"

namespace fem::material {

void Plasticity::update(const VolumeStrainLinearized&,
                        const Elasticity&,
                        const Precision*,
                        Precision*,
                        VolumeStressCauchy&,
                        Mat6&) const {
    logging::error(false, "Plasticity does not support linearized volume evaluation");
}

void Plasticity::recover(const VolumeStrainLinearized&,
                         const Elasticity&,
                         const Precision*,
                         VolumeStressCauchy&) const {
    logging::error(false, "Plasticity does not support linearized volume recovery");
}

} // namespace fem::material
