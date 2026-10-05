/**
 * @file dict.h
 * @brief Defines shared-object dictionaries with named or indexed storage.
 *
 * The model-data subsystem uses Dict to own shared definition objects and to
 * select an active object during parsing. String keys use an unordered map;
 * other keys address vector slots, including empty slots between sparse indices.
 * Creation, registration and activation have distinct effects on the active
 * pointer and are documented by the container contract.
 *
 * Stored objects own their domain behavior. String-keyed entries derive from
 * fem::Namable, while callers supply valid nonnegative indices for vector storage.
 * The complete template implementation is defined in this header.
 *
 * @see Dict
 * @see fem::Namable
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#include "../core/namable.h"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem {
namespace model {

namespace detail {

/**
 * @brief Defers the naming constraint while the stored object type is incomplete.
 *
 * This primary trait accepts forward declarations so model interfaces can declare
 * dictionaries before the referenced definition type is complete. The complete
 * specialization performs the actual inheritance check.
 *
 * @tparam T Type whose naming contract is queried.
 */
template<typename T, typename = void>
struct IsNamableOrIncomplete : std::true_type {};

/**
 * @brief Checks the naming base once the stored object type is complete.
 *
 * A successful sizeof(T) substitution selects the fem::Namable inheritance test.
 * The trait owns no state and supplies only the compile-time Boolean value used
 * by string-keyed Dict declarations.
 *
 * @tparam T Complete stored object type.
 */
template<typename T>
struct IsNamableOrIncomplete<T, std::void_t<decltype(sizeof(T))>>
    : std::bool_constant<std::is_base_of_v<fem::Namable, T>> {};

} // namespace detail

/**
 * @brief Owns shared objects addressed by names or nonnegative indices.
 *
 * String keys select unordered-map storage and require named objects. Other keys
 * select a vector and must be valid nonnegative indices; sparse creation grows the
 * vector with null slots. size() counts map keys or vector slots, not necessarily
 * nonnull objects. Iteration exposes the underlying representation, including
 * null entries and map key/pointer pairs.
 *
 * _cur retains the last activated or registered object. has_any() checks this
 * pointer, not whether storage is nonempty. Lookup and create() do not change it,
 * and remove() does not clear it, so an active object can outlive its registration.
 * Public storage and shared pointers allow mutation by callers; the dictionary
 * does not enforce domain-specific invariants or deep-copy the stored objects.
 *
 * activate() reuses an existing entry or creates and selects one. add() derives a
 * string key from the immutable object name and keeps an existing entry on a name
 * collision. Direct create() returns the new object even when string-key emplace
 * retains an old entry; callers should use activate() to obtain the registered
 * object. Indexed create() replaces the addressed slot. Shared ownership governs
 * lifetimes independently of registration.
 *
 * @tparam T Stored base object type.
 * @tparam Key std::string for named storage; otherwise a nonnegative vector index.
 */
template<typename T, typename Key = std::string>
struct Dict {
    // Shared ownership types for the dictionary and stored polymorphic objects.
    using Ptr  = std::shared_ptr<Dict<T, Key>>;
    using TPtr = std::shared_ptr<T>;

    // Owned registration storage: named map entries or indexed, possibly null slots.
    std::conditional_t<std::is_same_v<Key, std::string>, std::unordered_map<Key, TPtr>, std::vector<TPtr>> _data;

    static_assert(!std::is_same_v<Key, std::string> || detail::IsNamableOrIncomplete<T>::value,
                  "String-keyed Dict entries must derive from fem::Namable");

    // Shared active object; removal from storage does not reset this pointer.
    TPtr _cur = nullptr;

    // Named-key presence or nonnull indexed-slot presence; indices must be nonnegative.
    bool has(const Key& key) const {
        if constexpr (std::is_same_v<Key, std::string>) {
            return has_key(key);
        } else {
            return key < _data.size() && _data[key] != nullptr;
        }
    }

    // Check map-key presence only; indexed dictionaries always return false here.
    bool has_key(const Key& key) const {
        if constexpr (std::is_same_v<Key, std::string>) {
            return _data.find(key) != _data.end();
        }
        return false;
    }

    // Query/access the selected object, independently of whether storage is empty.
    bool has_any() const {
        return _cur != nullptr;
    }

    // Return the last object selected by activate() or named add().
    TPtr get() const {
        return _cur;
    }

    // Lookup without insertion or activation; missing entries return nullptr.
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

    // Reuse or create a registered object and retain it as the active selection.
    template<typename Derived = T, typename... Args>
    TPtr activate(const Key& key, Args... args) {
        static_assert(std::is_base_of_v<T, Derived>, "Derived must be derived from T");

        // Select only the requested object; get() itself does not update the active pointer.
        if (!has(key)) {
            _cur = create<Derived>(key, std::forward<Args>(args)...);
        } else {
            _cur = get(key);
        }
        return _cur;
    }

    // Register an already constructed named object. String-keyed dictionaries
    // derive the key from the object's immutable fem::Namable::name property so the
    // owning API does not need to accept the same name separately. A name collision
    // selects the existing registered object; null input leaves _cur unchanged.
    template<typename K = Key>
    std::enable_if_t<std::is_same_v<K, std::string>, TPtr> add(TPtr instance) {
        // Ignore null registration requests without disturbing the current selection.
        if (!instance) {
            return nullptr;
        }

        // Preserve the existing registration on collision and activate the stored pointer.
        const Key key = instance->name;
        auto [it, inserted] = _data.emplace(key, instance);
        (void) inserted;
        _cur = it->second;
        return _cur;
    }

    // Erase a map key or clear an indexed slot; preserve _cur and other shared owners.
    void remove(const Key& key) {
        if constexpr (std::is_same_v<Key, std::string>) {
            _data.erase(key);
        } else {
            if (key < _data.size()) {
                _data[key] = nullptr;
            }
        }
    }

    // Construct/register without activation. Named construction receives key first;
    // indexed construction receives only args and replaces an existing slot.
    template<typename Derived = T, typename... Args>
    TPtr create(const Key& key, Args... args) {
        static_assert(std::is_base_of_v<T, Derived>, "Derived must be derived from T");

        // Named constructors receive their identifier; indexed constructors do not.
        if constexpr (std::is_same_v<Key, std::string>) {
            auto instance = std::make_shared<Derived>(key, std::forward<Args>(args)...);
            _data.emplace(key, instance);
            return instance;
        } else {
            auto instance = std::make_shared<Derived>(std::forward<Args>(args)...);
            // Grow sparse indexed storage with null slots before assigning the object.
            if (key >= _data.size()) {
                _data.resize(key + 1);
            }
            _data[key] = instance;
            return instance;
        }
    }

    // Provides iterator support for range-based loops.
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
    // Count map keys or allocated vector slots, including any null entries.
    auto size() {
        return _data.size();
    }
    auto size() const {
        return _data.size();
    }
};
}    // namespace model
}    // namespace fem
