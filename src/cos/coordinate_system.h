/******************************************************************************
 * @file coordinate_system.hpp
 * @brief Defines the CoordinateSystem base class and RectangularSystem derived class for local
 *        coordinate systems used in rigid connectors. It provides methods for transforming
 *        points and vectors between global and local coordinate frames.
 *
 * @details The `CoordinateSystem` class serves as the base class for different local coordinate
 *          system implementations. It provides pure virtual methods for transforming points
 *          from global to local coordinates and vice versa. The `RectangularSystem` class
 *          defines a rectangular Cartesian coordinate system using two orthogonal vectors.
 *          This class is used for defining local reference frames in rigid connectors.
 *
 * @date Created on 08.10.2024
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include <Eigen/Dense>
#include <string>
#include <memory>

namespace fem {
namespace cos {

/**
 * @class CoordinateSystem
 * @brief Base class for representing local coordinate systems in FEM connectors.
 *
 * @details The CoordinateSystem class provides a virtual interface for defining
 *          local coordinate systems. Derived classes must implement the transformation
 *          functions for converting between global and local coordinate frames.
 */
class CoordinateSystem {
public:
    /// Default constructor
    CoordinateSystem() = default;

    /// Virtual destructor
    virtual ~CoordinateSystem() = default;

    /**
     * @brief Transform a point from global to local coordinates.
     * @param global_point Point in the global coordinate system.
     * @return Point in the local coordinate system.
     */
    virtual Eigen::Vector3d to_local(const Eigen::Vector3d& global_point) const = 0;

    /**
     * @brief Transform a point from local to global coordinates.
     * @param local_point Point in the local coordinate system.
     * @return Point in the global coordinate system.
     */
    virtual Eigen::Vector3d to_global(const Eigen::Vector3d& local_point) const = 0;

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
