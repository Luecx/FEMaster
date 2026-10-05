/**
 * @file castable.h
 * @brief Defines shared runtime casts for polymorphic FEMaster interfaces.
 *
 * The core Castable mixin provides mutable and const access to concrete types
 * through C++ runtime type information. Sections, elements, constitutive models
 * and load cases retain responsibility for their domain interfaces and state.
 *
 * Casts return non-owning pointers and preserve constness. A type mismatch yields
 * nullptr so callers can validate the required implementation before using it.
 *
 * @see Castable
 * @see Section
 * @see model::ElementInterface
 * @see material::Elasticity
 * @see loadcase::LoadCase
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

namespace fem {

/**
 * @brief Provides checked pointer access to polymorphic implementations.
 *
 * Public inheritance gives a domain interface both as<T>() overloads without
 * repeating dynamic_cast. The virtual destructor makes the mixin polymorphic
 * and permits destruction through its base pointer. Derived interfaces retain
 * all domain operations and storage; Castable owns no data or external objects.
 *
 * T must be a complete class type supported by dynamic_cast. Successful casts
 * point into the same object and do not change its state or lifetime. Failed
 * casts return nullptr. The const overload exposes only const access, and every
 * returned pointer remains valid only while the original object is alive.
 * Calling as<T>() requires a valid object; callers must check nullable source
 * pointers before invoking it.
 */
struct Castable {
    // Polymorphic lifetime management for every interface using this mixin.
    virtual ~Castable() = default;

    // Runtime access to mutable implementations; nullptr denotes a type mismatch.
    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }

    // Runtime access to const implementations with the same mismatch semantics.
    template<typename T>
    const T* as() const {
        return dynamic_cast<const T*>(this);
    }
};

} // namespace fem
