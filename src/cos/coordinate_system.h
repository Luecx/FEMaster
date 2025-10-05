/******************************************************************************
 * @file coordinate_system.h
 * @brief Declares the base interface for local coordinate systems.
 *
 * Coordinate systems provide transformations between global and local frames
 * used by constraints and elements that require orientation handling.
 *
 * @see src/cos/rectangular_system.h
 * @see src/cos/cylindrical_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include "../data/namable.h"

#include <memory>
#include <string>

namespace fem {
namespace cos {

using Basis = Mat3; ///< Convenience alias for a 3x3 basis matrix.

/******************************************************************************
 * @struct CoordinateSystem
 * @brief Abstract component that converts between global and local coordinates.
 ******************************************************************************/
struct CoordinateSystem : Namable {
    using Ptr = std::shared_ptr<CoordinateSystem>; ///< Shared pointer shorthand.

    explicit CoordinateSystem(const std::string& name = "") : Namable(name) {}
    virtual ~CoordinateSystem() = default;

    /******************************************************************************
     * @brief Converts a point from the global frame to the local frame.
     *
     * @param global_point Position expressed in global coordinates.
     * @return Vec3 Local coordinates.
     ******************************************************************************/
    virtual Vec3 to_local(const Vec3& global_point) const = 0;

    /******************************************************************************
     * @brief Converts a point from the local frame to the global frame.
     *
     * @param local_point Position expressed in local coordinates.
     * @return Vec3 Global coordinates.
     ******************************************************************************/
    virtual Vec3 to_global(const Vec3& local_point) const = 0;

    /******************************************************************************
     * @brief Returns the basis vectors of the local frame at a specific location.
     *
     * @param local_point Local coordinates where the basis is evaluated.
     * @return Basis Matrix containing the basis vectors as columns.
     ******************************************************************************/
    virtual Basis get_axes(const Vec3& local_point) const = 0;
};

} // namespace cos
} // namespace fem
