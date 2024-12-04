/******************************************************************************
 * @file collection.h
 * @brief Defines the Collection class for managing collections of generic items.
 *
 * @details The Collection class provides a templated interface for managing a collection
 *          of items, with optional sorting and duplicate management. It supports various
 *          operations such as adding items, accessing elements, and iterating over the
 *          collection. Special handling is included for non-sortable types to always allow
 *          duplicates and disable sorting.
 *
 * @date Created on 02.11.2024
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "namable.h"

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace fem::model {

/**
 * @class Collection
 * @brief Templated class for managing collections of items.
 *
 * @tparam T The type of items in the collection.
 *
 * @details The Collection class allows managing a list of items with optional support for
 *          sorting and duplicate removal. It provides an interface for adding items, iterating
 *          through them, and managing parent-child relationships between collections.
 */
template<typename T>
class Collection : public Namable{
public:
    using value_type = T; ///< The type of items stored in the collection.
    using Ptr = std::shared_ptr<Collection<T>>; ///< Shared pointer to a Collection.

    /**
     * @brief Constructor for the Collection.
     * @param p_name The name of the collection.
     * @param p_duplicates Whether duplicates are allowed.
     * @param p_sorted Whether the collection should maintain sorted order.
     */
    Collection(std::string p_name, bool p_duplicates = false, bool p_sorted = true)
        : Namable(p_name), _sorted(false), _duplicates(false) {
        sorted(p_sorted);
        duplicates(p_duplicates);
    }

    /**
     * @brief Enables or disables sorting for the collection.
     * @param p_sorted Whether to enable sorting.
     * @return Reference to the Collection for chaining.
     */
    Collection& sorted(bool p_sorted) {
        if constexpr (is_sortable) {
            if (p_sorted == _sorted) return *this;
            _sorted = p_sorted;
            if (_sorted) {
                std::sort(_data.begin(), _data.end());
            }
            if (!_sorted) {
                duplicates(true); // Ensure duplicates are allowed when sorting is disabled.
            }
        } else {
            _sorted = false; // Disable sorting for non-sortable types.
        }
        return *this;
    }

    /**
     * @brief Enables or disables duplicate removal for the collection.
     * @param p_duplicates Whether to allow duplicates.
     * @return Reference to the Collection for chaining.
     */
    Collection& duplicates(bool p_duplicates) {
        if constexpr (is_sortable) {
            if (p_duplicates == _duplicates) return *this;
            _duplicates = p_duplicates;
            if (!_duplicates) {
                sorted(true); // Sorting is required to remove duplicates.
                _data.erase(std::unique(_data.begin(), _data.end()), _data.end());
            }
        } else {
            _duplicates = true; // Always allow duplicates for non-sortable types.
        }
        return *this;
    }

    /**
     * @brief Sets the parent collection.
     * @param p_parent Pointer to the parent collection.
     */
    void set_parent(Ptr p_parent) {
        _parent = p_parent;
    }

    /**
     * @brief Adds an item to the collection.
     * @param item The item to add.
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
        if (_parent) _parent->add(item);
    }

    /**
     * @brief Adds all items from another collection.
     * @param items The collection of items to add.
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
        if (_parent) _parent->add(items);
    }

    // Iterators
    auto begin() { return _data.begin(); } ///< Iterator to the beginning of the collection.
    auto end() { return _data.end(); }   ///< Iterator to the end of the collection.
    auto begin() const { return _data.cbegin(); } ///< Const iterator to the beginning of the collection.
    auto end() const { return _data.cend(); }   ///< Const iterator to the end of the collection.


    /**
     * @brief Gets all items in the collection.
     * @return A const reference to the vector of items.
     */
    const std::vector<T>& data() const { return _data; }

    /**
     * @brief Gets the first item in the collection.
     * @return A const reference to the first item.
     */
    const T& first() const { return _data.front(); }

    /**
     * @brief Gets the last item in the collection.
     * @return A const reference to the last item.
     */
    const T& last() const { return _data.back(); }

    /**
     * @brief Gets the size of the collection.
     * @return The number of items in the collection.
     */
    size_t size() const { return _data.size(); }

    /**
     *
     */
    T& at(size_t index) {
        return _data[index];
    }

    T at(size_t index) const {
        return _data[index];
    }

    T operator[](size_t index) const {
        return _data[index];
    }

    T& operator[](size_t index) {
        return _data[index];
    }

    T operator()(size_t index) const {
        return _data[index];
    }

    T& operator()(size_t index) {
        return _data[index];
    }
protected:
    std::vector<T> _data; ///< The vector storing the items.
    Ptr _parent = nullptr; ///< Pointer to the parent collection (if any).

    bool _sorted; ///< Whether the collection maintains sorted order.
    bool _duplicates; ///< Whether duplicates are allowed.

private:
    /**
     * @brief Adds an item to a sorted collection.
     * @param item The item to add.
     */
    void add_sorted(const T& item) {
        auto it = std::lower_bound(_data.begin(), _data.end(), item);
        if (_duplicates || it == _data.end() || *it != item) {
            _data.insert(it, item);
        }
    }

    /**
     * @brief Adds an item to an unsorted collection.
     * @param item The item to add.
     */
    void add_unsorted(const T& item) {
        _data.push_back(item);
    }

    /**
     * @brief Adds multiple items to a sorted collection.
     * @param items The collection of items to add.
     */
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

    /**
     * @brief Adds multiple items to an unsorted collection.
     * @param items The collection of items to add.
     */
    void add_unsorted(const Collection<T>& items) {
        _data.insert(_data.end(), items._data.begin(), items._data.end());
    }

    // Trait to determine if a type is sortable
    static constexpr bool is_sortable = std::is_arithmetic<T>::value || std::is_pointer<T>::value;
};

}  // namespace fem::model

