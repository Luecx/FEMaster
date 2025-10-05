/**
 * @file namable.h
 * @brief Declares the lightweight naming mixin shared across model entities.
 *
 * The `Namable` struct stores a human-readable identifier that is used by
 * collections, material definitions, and other high-level FEM concepts. The
 * helper lives in the `fem::model` namespace but is re-exported into the global
 * namespace for backwards compatibility with existing includes.
 *
 * @see src/data/namable.cpp
 * @see src/data/collection.h
 * @author Finn Eggers
 * @date 07.03.2025
 */

#pragma once

#include <string>

namespace fem {
namespace model {

/**
 * @struct Namable
 * @brief Provides a shared immutable `name` property.
 */
struct Namable {
    const std::string name; ///< Human-readable identifier of the entity.

    /**
     * @brief Constructs the mixin with the supplied name.
     *
     * @param p_name Identifier assigned to the derived object.
     */
    explicit Namable(std::string p_name);
};

} // namespace model
} // namespace fem

using fem::model::Namable; ///< Retained alias for existing unqualified usage.

