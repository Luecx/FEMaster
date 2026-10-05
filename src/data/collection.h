/**
 * @file collection.h
 * @brief Defines named value collections with ordering and parent propagation.
 *
 * The model-data subsystem uses Collection as the storage and insertion-policy
 * base of typed regions. Values are owned in a vector; optional sorting and
 * duplicate suppression operate on arithmetic and raw pointer types. Insertions
 * can also propagate to a shared parent collection.
 *
 * Entity interpretation, region registration and aggregate-set selection belong
 * to Region and Sets. The template implementation is contained in this header.
 *
 * @see Collection
 * @see Region
 * @see Sets
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#include "../core/namable.h"

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace fem {
namespace model {

/**
 * @brief Owns named values and applies local insertion policies.
 *
 * Collection stores values by copy in a vector and owns its immutable name through
 * fem::Namable. Arithmetic and raw pointer types support ascending order and
 * duplicate suppression. Other types always use unsorted append operations and
 * allow duplicates; pointer ordering does not compare pointed-to objects.
 *
 * Disabling duplicates enables sorting and removes repeated values. Disabling
 * sorting enables duplicates. Mutable iterators and indexed references expose
 * storage directly, so callers must preserve these policy invariants themselves.
 * Indexed access is unchecked, and first()/last() require nonempty storage.
 * Insertions and policy changes may invalidate vector iterators and references.
 *
 * An optional shared parent receives every requested insertion, even if the local
 * collection suppresses it as a duplicate. Existing values are not propagated
 * when a parent is attached. Parent links must be acyclic; the class neither
 * checks cycles nor owns domain-specific entity objects beyond the stored T values.
 * Derived regions supply entity meaning and diagnostics without changing storage.
 *
 * @tparam T Stored value type; sorting is enabled only for arithmetic/raw pointers.
 */
template<typename T>
class Collection : public fem::Namable {
    public:
    // Stored value and shared collection ownership types.
    using value_type = T;
    using Ptr        = std::shared_ptr<Collection<T>>;

    // Construction with ordering and duplicate policies supported by T.
    Collection(std::string p_name, bool p_duplicates = false, bool p_sorted = true)
        : fem::Namable(std::move(p_name))
        , _sorted(false)
        , _duplicates(false) {
        // Establish ordering first, then enforce the requested duplicate policy.
        sorted(p_sorted);
        duplicates(p_duplicates);
    }

    // Enable ascending order, or allow duplicates when switching to unsorted storage.
    Collection& sorted(bool p_sorted) {
        // Select policy handling only when T supports the built-in ordering contract.
        if constexpr (is_sortable) {
            if (p_sorted == _sorted) {
                return *this;
            }
            _sorted = p_sorted;
            if (_sorted) {
                std::sort(_data.begin(), _data.end());
            }
            if (!_sorted) {
                duplicates(true);
            }
        } else {
            _sorted = false;
        }
        return *this;
    }

    // Allow duplicates, or sort and deduplicate the existing storage when disabling them.
    Collection& duplicates(bool p_duplicates) {
        // Select policy handling only when T supports the built-in ordering contract.
        if constexpr (is_sortable) {
            if (p_duplicates == _duplicates) {
                return *this;
            }
            _duplicates = p_duplicates;
            if (!_duplicates) {
                sorted(true);
                _data.erase(std::unique(_data.begin(), _data.end()), _data.end());
            }
        } else {
            _duplicates = true;
        }
        return *this;
    }

    // Attach a shared parent for future insertions; existing values are not copied.
    void set_parent(Ptr p_parent) {
        _parent = std::move(p_parent);
    }

    // Insert values according to local policies, then forward the request to the parent.
    void add(const T& item) {
        // Select policy handling only when T supports the built-in ordering contract.
        if constexpr (is_sortable) {
            if (_sorted) {
                add_sorted(item);
            } else {
                add_unsorted(item);
            }
        } else {
            add_unsorted(item);
        }
        // Forward the complete request even when the local insertion was suppressed.
        if (_parent) {
            _parent->add(item);
        }
    }

    // Insert another collection by value; source/target and parent storage must not alias.
    void add(const Collection<T>& items) {
        // Select policy handling only when T supports the built-in ordering contract.
        if constexpr (is_sortable) {
            if (_sorted) {
                add_sorted(items);
            } else {
                add_unsorted(items);
            }
        } else {
            add_unsorted(items);
        }
        // Forward the complete request even when the local insertion was suppressed.
        if (_parent) {
            _parent->add(items);
        }
    }

    // Iterator access to support range-based loops.
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

    // Returns the internal storage for read-only access.
    const std::vector<T>& data() const {
        return _data;
    }

    // Returns the first element. Undefined when the collection is empty.
    const T& first() const {
        return _data.front();
    }

    // Returns the last element. Undefined when the collection is empty.
    const T& last() const {
        return _data.back();
    }

    // Returns the current number of stored items.
    [[nodiscard]] size_t size() const {
        return _data.size();
    }

    // Unchecked mutable indexed access; callers must preserve ordering/uniqueness.
    T& at(size_t index) {
        return _data[index];
    }

    // Const indexed access returns a value copy; index must be less than size().
    T at(size_t index) const {
        return _data[index];
    }

    // Const subscript operator forwarding to `at`.
    T operator[](size_t index) const {
        return _data[index];
    }

    // Mutable subscript operator forwarding to `at`.
    T& operator[](size_t index) {
        return _data[index];
    }

    // Const call operator mirroring subscript semantics.
    T operator()(size_t index) const {
        return _data[index];
    }

    // Mutable call operator mirroring subscript semantics.
    T& operator()(size_t index) {
        return _data[index];
    }

    protected:
    // Owned value storage and shared parent receiving future insertion requests.
    std::vector<T> _data;
    Ptr            _parent = nullptr;

    // Effective insertion policies; uniqueness requires sorted storage.
    bool           _sorted;
    bool           _duplicates;

    private:
    // Inserts an item into the sorted storage while respecting duplicates.
    void add_sorted(const T& item) {
        // Find the ordered insertion position and reject an equal value if uniqueness is required.
        auto it = std::lower_bound(_data.begin(), _data.end(), item);
        if (_duplicates || it == _data.end() || *it != item) {
            _data.insert(it, item);
        }
    }

    // Appends an item without sorting.
    void add_unsorted(const T& item) {
        _data.push_back(item);
    }

    // Inserts items into a sorted collection.
    void add_sorted(const Collection<T>& items) {
        // Reuse the advancing lower-bound position only for an ordered source sequence.
        if (items._sorted) {
            auto it = _data.begin();
            for (const auto& item : items._data) {
                it = std::lower_bound(it, _data.end(), item);
                if (_duplicates || it == _data.end() || *it != item) {
                    it = _data.insert(it, item);
                    ++it;
                }
            }
        } else {
            for (const auto& item : items._data) {
                add_sorted(item);
            }
        }
    }

    // Appends items without sorting.
    void add_unsorted(const Collection<T>& items) {
        _data.insert(_data.end(), items._data.begin(), items._data.end());
    }

    // Sorting and duplicate comparisons are instantiated only for these value types.
    static constexpr bool is_sortable = std::is_arithmetic_v<T> || std::is_pointer_v<T>;
};
}    // namespace model
}    // namespace fem
