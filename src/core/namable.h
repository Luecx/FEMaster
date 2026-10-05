/**
 * @file namable.h
 * @brief Defines the header-only naming mixin shared by FEMaster entities.
 *
 * The core naming utility owns an immutable identifier for model definitions,
 * materials, profiles, coordinate systems and named collections. Dictionary
 * storage, name uniqueness and lookup remain responsibilities of the owning
 * containers; Namable neither registers objects nor validates their names.
 *
 * Construction moves the supplied string into the object. No separate
 * translation unit or dependency on model data is required.
 *
 * @see Namable
 * @see model::Dict
 * @see model::Collection
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#include <string>
#include <utility>

namespace fem {

/**
 * @brief Owns the immutable name of a persistent FEMaster definition.
 *
 * Derived entities reuse this core mixin to expose a common identifier without
 * coupling their domain data to dictionary or collection implementations. Each
 * object owns its string, and copying an object copies its name. The const name
 * cannot be reassigned, so assignment of the mixin is intentionally unavailable.
 *
 * Empty names are allowed. Naming rules and uniqueness within a definition scope
 * are enforced by callers and containers. The mixin has no runtime state or
 * virtual interface and is not intended for polymorphic deletion.
 */
struct Namable {
    // Owned identifier, fixed for the lifetime of this definition.
    const std::string name;

    // Take ownership of the supplied identifier without retaining caller storage.
    explicit Namable(std::string p_name) : name(std::move(p_name)) {}
};

} // namespace fem
