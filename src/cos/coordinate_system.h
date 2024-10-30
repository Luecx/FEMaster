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

#include "../core/types.h"
#include <Eigen/Dense>
#include <string>
#include <memory>

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
class CoordinateSystem {
public:

    Vec3 axis_1{0,0,0}; ///< First axis of the local coordinate system in global coordinates.
    Vec3 axis_2{0,0,0}; ///< Second axis of the local coordinate system in global coordinates.
    Vec3 axis_3{0,0,0}; ///< Third axis of the local coordinate system in global coordinates.

    /// Default constructor
    CoordinateSystem() = default;

    /// Constructor with three axes
    CoordinateSystem(Vec3 axis_1, Vec3 axis_2, Vec3 axis_3) {
        this->axis_1 = axis_1;
        this->axis_2 = axis_2;
        this->axis_3 = axis_3;
    };

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
     * @brief Get the name of the coordinate system.
     * @return The name of the coordinate system.
     */
    const std::string& get_name() const { return name_; }

    /**
     * @brief Set the name of the coordinate system.
     * @param name Name of the coordinate system.
     */
    void set_name(const std::string& name) { name_ = name; }

protected:
    std::string name_; ///< Name of the coordinate system.
};

using CoordinateSystemPtr = std::shared_ptr<CoordinateSystem>;

} // namespace fem
} // namespace fem::cos
