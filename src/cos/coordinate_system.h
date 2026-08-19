/**
 * @file coordinate_system.h
 * @brief Declares the base interface for local coordinate systems.
 *
 * Coordinate systems are definition objects and are deliberately independent
 * of the part/instance compilation lifecycle. They may therefore be registered
 * before or after `Model::compile()`. During compilation, an instance-local
 * section orientation is copied through `transformed()` so both directional
 * axes and spatial origins follow the rigid placement of the instance.
 *
 * @see src/cos/rectangular_system.h
 * @see src/cos/cylindrical_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/namable.h"

#include <memory>
#include <string>

namespace fem {
namespace cos {

using Basis = Mat3; ///< Convenience alias for a 3x3 basis matrix.

/**
 * @struct CoordinateSystem
 * @brief Abstract component that converts between global and local coordinates.
 *
 * Coordinate systems are immutable geometric definitions from the model point
 * of view. A section may reference one directly before compilation. When the
 * owning part is instantiated, `Model::compile()` requests a rigidly transformed
 * copy instead of mutating the shared source definition.
 */
struct CoordinateSystem : Namable {
    using Ptr = std::shared_ptr<CoordinateSystem>; ///< Shared pointer shorthand.

    explicit CoordinateSystem(const std::string& name = "") : Namable(name) {}
    virtual ~CoordinateSystem() = default;

    /**
     * @brief Converts a point from the global frame to the local frame.
     *
     * @param global_point Position expressed in global coordinates.
     * @return Vec3 Local coordinates.
     */
    virtual Vec3 to_local(const Vec3& global_point) const = 0;

    /**
     * @brief Converts a point from the local frame to the global frame.
     *
     * @param local_point Position expressed in local coordinates.
     * @return Vec3 Global coordinates.
     */
    virtual Vec3 to_global(const Vec3& local_point) const = 0;

    /**
     * @brief Returns the basis vectors of the local frame at a specific location.
     *
     * @param local_point Local coordinates where the basis is evaluated.
     * @return Basis Matrix containing the basis vectors as columns.
     */
    virtual Basis get_axes(const Vec3& local_point) const = 0;

    /**
     * @brief Creates a copy embedded by one rigid instance transformation.
     *
     * The transformation follows the same convention as `Instance`: a source
     * point `x` becomes `rotation * x + translation`. Direction vectors are
     * affected only by `rotation`; coordinate systems with a spatial origin
     * must transform that origin by both rotation and translation.
     *
     * @param rotation Proper orthonormal instance rotation.
     * @param translation Instance translation in global coordinates.
     * @return Independent coordinate-system definition in the instance frame.
     */
    virtual Ptr transformed(const Mat3& rotation, const Vec3& translation) const = 0;
};
} // namespace cos
} // namespace fem
