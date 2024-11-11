//
// Created by Finn Eggers on 02.11.24.
//

#ifndef FEMASTER_COLLECTION_H
#define FEMASTER_COLLECTION_H

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <memory>

namespace fem::model {

template<typename T>
struct Collection {
    using value_type = T;
    using Ptr = std::shared_ptr<Collection<T>>;

    protected:
    std::vector<T> _data;
    Ptr _parent = nullptr;

    const std::string _name;

    bool _sorted;
    bool _duplicates;

    public:
    Collection(std::string p_name, bool p_duplicates = false, bool p_sorted = true)
        : _name(std::move(p_name)) {
        sorted(p_sorted);
        duplicates(p_duplicates);
    }

    Collection<T>& sorted(bool p_sorted) {
        if (p_sorted == _sorted) return *this;
        _sorted = p_sorted;
        if (_sorted) {
            std::sort(_data.begin(), _data.end());
        }
        if (!_sorted) {
            duplicates(true);
        }
        return *this;
    }

    Collection<T>& duplicates(bool p_duplicates) {
        if (p_duplicates == _duplicates) return *this;
        _duplicates = p_duplicates;
        if (!_duplicates) {
            sorted(true);
            _data.erase(std::unique(_data.begin(), _data.end()), _data.end());
        }
        return *this;
    }

    void set_parent(Ptr p_parent) {
        _parent = p_parent;
    }

    void add(const T& item) {
        if (_sorted) {
            add_sorted(item);
        } else {
            add_unsorted(item);
        }
        if (_parent) _parent->add(item);
    }

    void add(const Collection<T>& items) {
        if (_sorted) {
            add_sorted(items);
        } else {
            add_unsorted(items);
        }
        if (_parent) _parent->add(items);
    }

    // iterators
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }
    auto begin() const { return _data.cbegin(); }
    auto end() const { return _data.cend(); }

    // accessors
    const std::string& name() const { return _name; }
    const std::vector<T>& data() const { return _data; }

    // getting the first element
    const T& first() const { return _data.front(); }

    // getting the last element
    const T& last() const { return _data.back(); }

    // size
    size_t size() const { return _data.size(); }

    private:
    // Helper function for sorted insertion of a single item
    void add_sorted(const T& item) {
        auto it = std::lower_bound(_data.begin(), _data.end(), item);
        if (_duplicates || it == _data.end() || *it != item) {
            _data.insert(it, item);
        }
    }

    // Helper function for unsorted insertion of a single item
    // Appends the item at the end, assuming unsorted and duplicates allowed
    void add_unsorted(const T& item) {
        _data.push_back(item);
    }

    // Helper function for sorted insertion of multiple items from another Collection
    void add_sorted(const Collection<T>& items) {
        if (items._sorted) {
            auto it = _data.begin();
            for (const auto& item : items._data) {
                it = std::lower_bound(it, _data.end(), item);
                if (_duplicates || it == _data.end() || *it != item) {
                    it = _data.insert(it, item);
                    ++it;  // Move iterator past the newly inserted item
                }
            }
        } else {
            for (const auto& item : items._data) {
                add_sorted(item);
            }
        }
    }

    // Helper function for unsorted insertion of multiple items from another Collection
    void add_unsorted(const Collection<T>& items) {
        _data.insert(_data.end(), items._data.begin(), items._data.end());
    }
};

}
#endif    // FEMASTER_COLLECTION_H
