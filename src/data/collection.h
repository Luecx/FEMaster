/**
 * @file collection.h
 * @brief Declares a generic container that manages named collections of items.
 *
 * The `Collection` template supplements `std::vector` with optional sorting,
 * duplicate handling, and parent propagation. It forms the basis for higher
 * level FEM sets such as node or element regions.
 *
 * @see src/data/collection.cpp
 * @see src/data/sets.h
 */

#pragma once

#include "namable.h"

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace fem {
namespace model {

/**
 * @class Collection
 * @brief Manages a list of items with optional sorting and duplicate control.
 *
 * @tparam T Value type stored inside the collection.
 */
template<typename T>
class Collection : public Namable {
public:
    using value_type = T;                      ///< Value stored in the collection.
    using Ptr = std::shared_ptr<Collection<T>>; ///< Shared pointer alias for derived collections.

    /**
     * @brief Constructs a collection with the supplied name and policies.
     *
     * Sorting and duplicate removal are only available when `T` is sortable
     * (either arithmetic or pointer types).
     */
    Collection(std::string p_name, bool p_duplicates = false, bool p_sorted = true)
        : Namable(std::move(p_name)), _sorted(false), _duplicates(false) {
        sorted(p_sorted);
        duplicates(p_duplicates);
    }

    /**
     * @brief Enables or disables automatic sorting of inserted items.
     */
    Collection& sorted(bool p_sorted) {
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

    /**
     * @brief Configures whether the collection suppresses duplicates.
     */
    Collection& duplicates(bool p_duplicates) {
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

    /// Assigns a parent collection that mirrors insertions.
    void set_parent(Ptr p_parent) { _parent = std::move(p_parent); }

    /**
     * @brief Inserts a single item respecting the configured policies.
     */
    void add(const T& item) {
        if constexpr (is_sortable) {
            if (_sorted) {
                add_sorted(item);
            } else {
                add_unsorted(item);
            }
        } else {
            add_unsorted(item);
        }
        if (_parent) {
            _parent->add(item);
        }
    }

    /**
     * @brief Inserts all items from another collection.
     */
    void add(const Collection<T>& items) {
        if constexpr (is_sortable) {
            if (_sorted) {
                add_sorted(items);
            } else {
                add_unsorted(items);
            }
        } else {
            add_unsorted(items);
        }
        if (_parent) {
            _parent->add(items);
        }
    }

    /// Iterator access to support range-based loops.
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }
    auto begin() const { return _data.cbegin(); }
    auto end() const { return _data.cend(); }

    /// Returns the internal storage for read-only access.
    const std::vector<T>& data() const { return _data; }

    /// Returns the first element. Undefined when the collection is empty.
    const T& first() const { return _data.front(); }

    /// Returns the last element. Undefined when the collection is empty.
    const T& last() const { return _data.back(); }

    /// Returns the current number of stored items.
    [[nodiscard]] size_t size() const { return _data.size(); }

    /// Mutable indexed access.
    T& at(size_t index) { return _data[index]; }

    /// Const indexed access.
    T at(size_t index) const { return _data[index]; }

    /// Const subscript operator forwarding to `at`.
    T operator[](size_t index) const { return _data[index]; }

    /// Mutable subscript operator forwarding to `at`.
    T& operator[](size_t index) { return _data[index]; }

    /// Const call operator mirroring subscript semantics.
    T operator()(size_t index) const { return _data[index]; }

    /// Mutable call operator mirroring subscript semantics.
    T& operator()(size_t index) { return _data[index]; }

protected:
    std::vector<T> _data; ///< Backing container storing the items.
    Ptr _parent = nullptr; ///< Optional parent that mirrors insertions.

    bool _sorted;     ///< Indicates whether the collection maintains sorted order.
    bool _duplicates; ///< Indicates whether duplicates are allowed.

private:
    /// Inserts an item into the sorted storage while respecting duplicates.
    void add_sorted(const T& item) {
        auto it = std::lower_bound(_data.begin(), _data.end(), item);
        if (_duplicates || it == _data.end() || *it != item) {
            _data.insert(it, item);
        }
    }

    /// Appends an item without sorting.
    void add_unsorted(const T& item) { _data.push_back(item); }

    /// Inserts items into a sorted collection.
    void add_sorted(const Collection<T>& items) {
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

    /// Appends items without sorting.
    void add_unsorted(const Collection<T>& items) {
        _data.insert(_data.end(), items._data.begin(), items._data.end());
    }

    /// Trait indicating whether sorting operations are available for `T`.
    static constexpr bool is_sortable = std::is_arithmetic_v<T> || std::is_pointer_v<T>;
};

} // namespace model
} // namespace fem

