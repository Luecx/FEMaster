/**
 * @file coordinate_system.h
 * @brief Defines the common interface for local geometric coordinate systems.
 *
 * The coordinate-system subsystem supplies point mappings and local orthonormal
 * bases for material orientations and constraint transformations. Basis columns
 * are local unit directions expressed in global Cartesian coordinates. Sections
 * and constraints apply these bases to their own vectors, tensors and DOFs.
 *
 * Coordinate systems are persistent definitions independent of model compilation.
 * An instance placement creates an independent transformed definition instead of
 * modifying a shared source. Concrete systems determine whether a spatial origin
 * or only an orientation belongs to their definition.
 *
 * @see CoordinateSystem
 * @see RectangularSystem
 * @see CylindricalSystem
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_eig.h"
#include "../core/namable.h"

#include <memory>
#include <string>

namespace fem {
namespace cos {

// Local unit directions expressed in global Cartesian coordinates as columns.
using Basis = Mat3;

/**
 * @brief Named polymorphic interface for point mappings and local bases.
 *
 * The base owns only its immutable name through fem::Namable. Derived systems own
 * their geometric definition and implement point mappings, basis evaluation and
 * copying under rigid placement. Definitions carry no constitutive history or
 * solver state and are shared through Ptr without mutation during evaluation.
 *
 * For B = get_axes(local_point), columns are the local unit directions in global
 * Cartesian coordinates. Thus B maps local vector components to global components
 * and B.transpose() maps global vector components to local components. Point
 * mappings additionally depend on the concrete coordinate representation and its
 * origin; cylindrical point coordinates are not vector components in this basis.
 *
 * transformed() follows the instance convention x' = rotation * x + translation.
 * The supplied rotation must be proper orthonormal. Direction vectors rotate,
 * while any spatial origin also translates. The copy retains the semantic name;
 * the shared source definition remains unchanged.
 */
struct CoordinateSystem : fem::Namable {
    // Shared ownership of persistent coordinate-system definitions referenced by
    // sections and other model entities. Multiple consumers may reuse the same
    // definition; const evaluation does not modify its geometric data.
    using Ptr = std::shared_ptr<CoordinateSystem>;

    // Construction and polymorphic destruction.
    explicit CoordinateSystem(const std::string& name = "") : fem::Namable(name) {}
    virtual ~CoordinateSystem() = default;

    // Point mappings between global Cartesian and concrete local coordinates.
    // to_local() expresses a global point in the local representation; to_global()
    // reconstructs a global point and includes any spatial origin. Rectangular
    // systems use Cartesian components, cylindrical systems their documented
    // radial/angular/axial representation with angles in radians. Whether these
    // operations are inverses depends on the concrete coordinate conventions.
    // Vector and tensor components are transformed using get_axes() instead.
    virtual Vec3 to_local(const Vec3& global_point) const = 0;
    virtual Vec3 to_global(const Vec3& local_point) const = 0;

    // Evaluate the orthonormal local basis at a point in local coordinates.
    // The returned columns are local unit directions expressed in global Cartesian
    // coordinates: B maps local vector components to global components, and
    // B.transpose() performs the reverse mapping. Rectangular bases are constant;
    // cylindrical radial/tangential directions depend on the supplied azimuth.
    virtual Basis get_axes(const Vec3& local_point) const = 0;

    // Create an independent definition under the rigid instance placement
    // x' = rotation * x + translation, retaining the name and leaving this object
    // unchanged. rotation must be proper orthonormal. Unit directions rotate;
    // a spatial origin also translates. Orientation-only systems ignore translation.
    // The returned shared pointer owns the newly created concrete definition.
    virtual Ptr transformed(const Mat3& rotation, const Vec3& translation) const = 0;
};
} // namespace cos
} // namespace fem
