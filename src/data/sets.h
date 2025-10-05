/******************************************************************************
 * @file sets.h
 * @brief Declares the `Sets` helper for managing named collections of regions.
 *
 * A `Sets` instance acts like a registry of named `Collection` objects. It
 * keeps track of an optional "all" collection that aggregates every inserted
 * value and provides convenience helpers for activation and iteration.
 *
 * @see src/data/sets.cpp
 * @see src/data/collection.h
 ******************************************************************************/

#pragma once

#include "collection.h"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace fem {
namespace model {

#define SET_NODE_ALL "NALL" ///< Convenience literal describing the global node set.
#define SET_ELEM_ALL "EALL" ///< Convenience literal describing the global element set.
#define SET_SURF_ALL "SFALL" ///< Convenience literal describing the global surface set.

/******************************************************************************
 * @struct Sets
 * @brief Manages a family of named collections and an optional aggregate set.
 *
 * @tparam T Collection type derived from `Collection<ValueType>`.
 ******************************************************************************/
template<typename T>
struct Sets {
    using ValueType = typename T::value_type; ///< Values stored inside each collection.
    using TPtr = typename T::Ptr;             ///< Shared pointer alias for collections.
    using Key = std::string;                  ///< Key type used to address collections.

    static_assert(std::is_base_of_v<Collection<ValueType>, T>,
                  "T must derive from Collection<ValueType>");

    Key _all_key; ///< Name of the aggregate collection.

    TPtr _all; ///< Aggregate collection storing every inserted value.
    TPtr _cur; ///< Currently active collection targeted by `add` operations.
    std::unordered_map<Key, TPtr> _data; ///< Map of named collections.

    /// Returns `true` if `name` refers to an existing collection.
    bool has(const Key& name) { return has_key(name); }

    /// Checks existence purely via lookup. Identical to `has` for string keys.
    bool has_key(const Key& name) { return _data.find(name) != _data.end(); }

    /// Returns whether the aggregate collection exists.
    bool has_all() { return _all != nullptr; }

    /// Returns whether any collection is currently active.
    bool has_any() { return _cur != nullptr; }

    /// Provides access to the aggregate collection.
    TPtr all() { return _all; }

    /// Returns the currently active collection.
    TPtr get() { return _cur; }

    /// Retrieves the collection associated with `name` (may return `nullptr`).
    TPtr get(const Key& name) { return _data[name]; }

    /******************************************************************************
     * @brief Constructs the registry and optionally creates the aggregate set.
     ******************************************************************************/
    explicit Sets(const Key& all_key = "")
        : _all_key(all_key) {
        if (!_all_key.empty()) {
            _all = create(_all_key);
            _all->sorted(true);
            _all->duplicates(false);
        }
    }

    /******************************************************************************
     * @brief Activates the collection associated with `name`, creating it when needed.
     *
     * Passing an empty `name` activates the aggregate collection when available.
     ******************************************************************************/
    template<typename... Args>
    TPtr activate(const Key& name, Args... c) {
        Key key = name;
        if (key.empty()) {
            key = _all_key;
        }
        if (key == _all_key) {
            _cur = _all;
        } else {
            if (!has_key(key)) {
                _cur = create(name, std::forward<Args>(c)...);
            } else {
                _cur = _data[name];
            }
        }
        return _cur;
    }

    /******************************************************************************
     * @brief Adds a single value to the active and aggregate collections.
     ******************************************************************************/
    void add(const ValueType& item) {
        if (_cur) {
            _cur->add(item);
        }
        if (_all) {
            _all->add(item);
        }
    }

    /******************************************************************************
     * @brief Adds an arithmetic range [first, last] with the provided `step`.
     ******************************************************************************/
    template<typename U = ValueType>
    std::enable_if_t<std::is_integral_v<U>> add(U first, U last, U step) {
        for (U value = first; value <= last; value += step) {
            add(value);
        }
    }

    /// Iterator access to the underlying associative container.
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }

private:
    /******************************************************************************
     * @brief Creates a new collection and registers it under `name`.
     ******************************************************************************/
    template<typename... Args>
    TPtr create(const Key& name, Args... c) {
        auto collection = std::make_shared<T>(name, std::forward<Args>(c)...);
        if (_all) {
            collection->set_parent(_all);
        }
        _data.emplace(name, collection);
        return collection;
    }
};

} // namespace model
} // namespace fem

