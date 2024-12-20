/******************************************************************************
 * @file coordinate_system.hpp
 * @brief Defines the CoordinateSystem base class and RectangularSystem derived class for local
 *        coordinate systems used in rigid _connectors. It provides methods for transforming
 *        points and vectors between global and local coordinate frames.
 *
 * @details The `CoordinateSystem` class serves as the base class for different local coordinate
 *          system implementations. It provides pure virtual methods for transforming points
 *          from global to local coordinates and vice versa. The `RectangularSystem` class
 *          defines a rectangular Cartesian coordinate system using two orthogonal vectors.
 *          This class is used for defining local reference frames in rigid _connectors.
 *
 * @date Created on 08.10.2024
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include "../data/namable.h"

#include <Eigen/Dense>
#include <memory>
#include <string>

namespace fem {
namespace cos {

/**
 * @class CoordinateSystem
 * @brief Base class for representing local coordinate systems in FEM _connectors.
 *
 * @details The CoordinateSystem class provides a virtual interface for defining
 *          local coordinate systems. Derived classes must implement the transformation
 *          functions for converting between global and local coordinate frames.
 */
struct CoordinateSystem : Namable {
public:
    using Ptr = std::shared_ptr<CoordinateSystem>;

    /// Default constructor
    CoordinateSystem(const std::string& name="") : Namable(name) {};

    /// Virtual destructor
    virtual ~CoordinateSystem() = default;

    /**
     * @brief Transform a point from global to local coordinates.
     * @param global_point Point in the global coordinate system.
     * @return Point in the local coordinate system.
     */
    virtual Vec3 to_local(const Vec3& global_point) const = 0;

    /**
     * @brief Transform a point from local to global coordinates.
     * @param local_point Point in the local coordinate system.
     * @return Point in the global coordinate system.
     */
    virtual Vec3 to_global(const Vec3& local_point) const = 0;

    /**
     * @brief gets the axes at a specified local point
     */
    virtual StaticMatrix<3,3> get_axes(const Vec3& local_point) const = 0;
};


} // namespace fem
} // namespace fem::cos
