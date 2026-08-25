/**
 * @file section_point_mass.h
 * @brief Declares concentrated point mass, inertia and ground impedance properties.
 *
 * `PointMassSection` is the physical property assigned to zero-dimensional
 * `model::PointElement` topology. The element owns identity and ELSET membership;
 * the section owns concentrated translational mass, rotary inertia and optional
 * translational/rotational stiffness and viscous damping against ground.
 *
 * The section intentionally has no material. Its values are assembled directly
 * by the point element and therefore remain independent of continuum density or
 * constitutive laws.
 *
 * @see src/model/element/point.h
 * @see src/section/section.h
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "section.h"

namespace fem {

/**
 * @struct PointMassSection
 * @brief Concentrated mass, rotary inertia and diagonal ground impedance.
 *
 * One section may be assigned to any number of one-node point elements through
 * its element region. Translational mass is isotropic, rotary inertia is stored
 * per rotational axis, and stiffness/damping vectors act independently against
 * ground in the six generalized nodal directions.
 */
struct PointMassSection : Section {
    using Ptr = std::shared_ptr<PointMassSection>;

    Precision mass_                     = Precision(0); ///< Translational point mass on ux, uy and uz.
    Vec3      rotary_inertia_           = Vec3::Zero(); ///< Rotary inertia on rx, ry and rz.
    Vec3      spring_constants_         = Vec3::Zero(); ///< Translational ground stiffness on ux, uy and uz.
    Vec3      rotary_spring_constants_  = Vec3::Zero(); ///< Rotational ground stiffness on rx, ry and rz.
    Vec3      damping_constants_        = Vec3::Zero(); ///< Translational viscous damping against ground.
    Vec3      rotary_damping_constants_ = Vec3::Zero(); ///< Rotational viscous damping against ground.

    PointMassSection(model::ElementRegion::Ptr region,
                     Precision                 mass,
                     Vec3                      rotary_inertia = Vec3::Zero(),
                     Vec3                      spring_constants = Vec3::Zero(),
                     Vec3                      rotary_spring_constants = Vec3::Zero(),
                     Vec3                      damping_constants = Vec3::Zero(),
                     Vec3                      rotary_damping_constants = Vec3::Zero());

    void info() override;
    std::string str() const override;
};

} // namespace fem
