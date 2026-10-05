/**
 * @file sets.h
 * @brief Defines named collection registries and optional aggregate sets.
 *
 * The model-data subsystem uses Sets to select a named region during input and
 * propagate inserted entity values to an optional all-entities collection. Named
 * collections are shared objects; their ordering and duplicate policies belong
 * to Collection. Region supplies the entity kind rather than this registry.
 *
 * The complete registry implementation is defined in this header. It manages
 * activation and parent links, while callers own model topology and identifier
 * validity.
 *
 * @see Sets
 * @see Collection
 * @see Region
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#include "collection.h"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace fem {
namespace model {

// Conventional aggregate name for node collections.
#define SET_NODE_ALL "NALL"
// Conventional aggregate name for element collections.
#define SET_ELEM_ALL "EALL"
// Conventional aggregate name for surface collections.
#define SET_SURF_ALL "SFALL"
// Conventional aggregate name for line collections.
#define SET_LINE_ALL "LALL"

/**
 * @brief Owns a registry of named collections with an active and aggregate set.
 *
 * _data holds shared collections keyed by name. _cur selects the insertion target,
 * and _all optionally aggregates requests across targets. Construction with a
 * nonempty aggregate name creates that collection and requests sorted unique
 * storage; the effective policy still depends on the collection value type.
 * New named collections receive _all as a shared parent when it exists.
 *
 * activate() maps an empty name to the aggregate name and creates missing named
 * collections. With no aggregate, empty-name activation selects nullptr. add()
 * updates the active collection and independently the aggregate. Parent forwarding
 * can therefore repeat the aggregate insertion request; sorted unique aggregate
 * storage suppresses repeated arithmetic/pointer values. Existing values are not
 * retroactively propagated when a parent is attached.
 *
 * get(name) uses map subscript access: a missing name inserts a null entry without
 * creating a collection or changing _cur. has() reports key presence even for such
 * an entry, and activation of it selects nullptr. Public pointers and storage may
 * be changed by callers; this registry does not validate entity identifiers.
 * Iteration returns unordered map entries. Parent ownership must remain acyclic.
 *
 * @tparam T Named collection type deriving from Collection<T::value_type>.
 */
template<typename T>
struct Sets {
    // Collection values, shared ownership and named lookup types.
    using ValueType = typename T::value_type;
    using TPtr = typename T::Ptr;
    using Key = std::string;

    static_assert(std::is_base_of_v<Collection<ValueType>, T>,
                  "T must derive from Collection<ValueType>");

    // Configured aggregate name; an empty string disables aggregate creation.
    Key _all_key;

    // Shared aggregate/active collections and owned named registrations.
    TPtr _all;
    TPtr _cur;
    std::unordered_map<Key, TPtr> _data;

    // Named-key presence queries; null map entries still count as registered keys.
    bool has(const Key& name) { return has_key(name); }

    // Checks existence purely via lookup. Identical to `has` for string keys.
    bool has_key(const Key& name) { return _data.find(name) != _data.end(); }

    // Returns whether the aggregate collection exists.
    bool has_all() { return _all != nullptr; }
    bool has_all() const { return _all != nullptr; }

    // Returns whether any collection is currently active.
    bool has_any() { return _cur != nullptr; }

    // Provides access to the aggregate collection.
    TPtr all() { return _all; }
    TPtr all() const { return _all; }

    // Returns the currently active collection.
    TPtr get() { return _cur; }

    // Lookup through map subscript; missing names insert a null entry without activation.
    TPtr get(const Key& name) { return _data[name]; }

    // Construction with an optional aggregate name.
    explicit Sets(const Key& all_key = "")
        : _all_key(all_key) {
        // Request sorted unique aggregate storage when an aggregate name is configured.
        if (!_all_key.empty()) {
            _all = create(_all_key);
            _all->sorted(true);
            _all->duplicates(false);
        }
    }

    // Select/create a named collection; an empty name selects the optional aggregate.
    template<typename... Args>
    TPtr activate(const Key& name, Args... c) {
        // Resolve the empty-name shorthand before selecting aggregate or named storage.
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

    // Forward one insertion request to both active and aggregate collections.
    void add(const ValueType& item) {
        // Insert into the selected collection; its parent link may already update _all.
        if (_cur) {
            _cur->add(item);
        }
        // Independently maintain the aggregate, including when no collection is active.
        if (_all) {
            _all->add(item);
        }
    }

    // Insert an inclusive ascending integral range. step must be positive and
    // all increments representable; the implementation does not check these conditions.
    template<typename U = ValueType>
    std::enable_if_t<std::is_integral_v<U>> add(U first, U last, U step) {
        // Include the upper bound when reached by positive, representable increments.
        for (U value = first; value <= last; value += step) {
            add(value);
        }
    }

    // Iterator access to the underlying associative container.
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }
    auto begin() const { return _data.cbegin(); }
    auto end() const { return _data.cend(); }

private:
    // Create a named collection and attach the existing aggregate as its shared parent.
    template<typename... Args>
    TPtr create(const Key& name, Args... c) {
        // Own the new collection and connect future insertions to the aggregate.
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
