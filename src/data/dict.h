/**
 * @file dict.h
 * @brief Declares a heterogeneous dictionary for FEM entities.
 *
 * The `Dict` template offers a lightweight associative container that stores
 * shared pointers either in an unordered map (string keys) or in a vector
 * (integral keys). It is primarily used to manage lazily-created model
 * components such as regions and fields.
 *
 * @see src/data/dict.cpp
 * @see src/data/collection.h
 */

#pragma once

#include "namable.h"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem {
namespace model {

namespace detail {

template<typename T, typename = void>
struct IsNamableOrIncomplete : std::true_type {};

template<typename T>
struct IsNamableOrIncomplete<T, std::void_t<decltype(sizeof(T))>>
    : std::bool_constant<std::is_base_of_v<Namable, T>> {};

} // namespace detail

/**
 * @struct Dict
 * @brief Stores shared pointers keyed either by string or by index.
 *
 * @tparam T   Stored type.
 * @tparam Key Key type (`std::string` for associative access, otherwise assumed
 *             to be an indexable integral type).
 */
template<typename T, typename Key = std::string>
struct Dict {
    using Ptr  = std::shared_ptr<Dict<T, Key>>;    ///< Shared pointer alias used by callers.
    using TPtr = std::shared_ptr<T>;               ///< Shared pointer alias for stored values.

    /// Backing container, selected depending on the key type.
    std::conditional_t<std::is_same_v<Key, std::string>, std::unordered_map<Key, TPtr>, std::vector<TPtr>> _data;

    static_assert(!std::is_same_v<Key, std::string> || detail::IsNamableOrIncomplete<T>::value,
                  "String-keyed Dict entries must derive from Namable");

    TPtr _cur = nullptr;    ///< Tracks the most recently accessed entry.

    /**
     * @brief Checks whether an entry exists for `key`.
     */
    bool has(const Key& key) const {
        if constexpr (std::is_same_v<Key, std::string>) {
            return has_key(key);
        } else {
            return key < _data.size() && _data[key] != nullptr;
        }
    }

    /**
     * @brief Checks existence using only key look-up (string specialisation).
     */
    bool has_key(const Key& key) const {
        if constexpr (std::is_same_v<Key, std::string>) {
            return _data.find(key) != _data.end();
        }
        return false;
    }

    /// Returns whether the dictionary currently holds any entry.
    bool has_any() const {
        return _cur != nullptr;
    }

    /// Returns the last activated/created entry.
    TPtr get() const {
        return _cur;
    }

    /**
     * @brief Retrieves the entry associated with `key`.
     *
     * @return Shared pointer or `nullptr` when the entry is missing.
     */
    TPtr get(const Key& key) const {
        if (!has(key)) {
            return nullptr;
        }

        if constexpr (std::is_same_v<Key, std::string>) {
            return _data.at(key);
        } else {
            return key < _data.size() ? _data[key] : nullptr;
        }
    }

    /**
     * @brief Activates the entry at `key`, creating it if necessary.
     *
     * @tparam Derived Concrete type deriving from `T`.
     * @tparam Args    Constructor argument types.
     */
    template<typename Derived = T, typename... Args>
    TPtr activate(const Key& key, Args... args) {
        static_assert(std::is_base_of_v<T, Derived>, "Derived must be derived from T");

        if (!has(key)) {
            _cur = create<Derived>(key, std::forward<Args>(args)...);
        } else {
            _cur = get(key);
        }
        return _cur;
    }

    // Register an already constructed named object. String-keyed dictionaries
    // derive the key from the object's immutable Namable::name property so the
    // owning API does not need to accept the same name separately.
    template<typename K = Key>
    std::enable_if_t<std::is_same_v<K, std::string>, TPtr> add(TPtr instance) {
        if (!instance) {
            return nullptr;
        }

        const Key key = instance->name;
        auto [it, inserted] = _data.emplace(key, instance);
        (void) inserted;
        _cur = it->second;
        return _cur;
    }

    /**
     * @brief Removes the entry identified by `key`.
     */
    void remove(const Key& key) {
        if constexpr (std::is_same_v<Key, std::string>) {
            _data.erase(key);
        } else {
            if (key < _data.size()) {
                _data[key] = nullptr;
            }
        }
    }

    /**
     * @brief Creates and stores a new instance associated with `key`.
     *
     * When `Key` is a string, the constructor forwards the key to the stored
     * type if it inherits from `Namable`.
     */
    template<typename Derived = T, typename... Args>
    TPtr create(const Key& key, Args... args) {
        static_assert(std::is_base_of_v<T, Derived>, "Derived must be derived from T");

        if constexpr (std::is_same_v<Key, std::string>) {
            auto instance = std::make_shared<Derived>(key, std::forward<Args>(args)...);
            _data.emplace(key, instance);
            return instance;
        } else {
            auto instance = std::make_shared<Derived>(std::forward<Args>(args)...);
            if (key >= _data.size()) {
                _data.resize(key + 1);
            }
            _data[key] = instance;
            return instance;
        }
    }

    /// Provides iterator support for range-based loops.
    auto begin() {
        return _data.begin();
    }
    auto end() {
        return _data.end();
    }
    auto begin() const {
        return _data.cbegin();
    }
    auto end() const {
        return _data.cend();
    }
    auto size() {
        return _data.size();
    }
    auto size() const {
        return _data.size();
    }
};
}    // namespace model
}    // namespace fem
