//
// Created by Finn Eggers on 02.11.24.
//

#ifndef FEMASTER_SETS_H
#define FEMASTER_SETS_H

#include <type_traits>
#include <memory>
#include "collection.h"
#include <unordered_map>

namespace fem::model {

#define SET_NODE_ALL "NALL"
#define SET_ELEM_ALL "EALL"
#define SET_SURF_ALL "SFALL"

template <typename T>
struct Sets {

    // Deduce the contained type by extracting it from Collection<T>::_data
    using ValueType = typename T::value_type;
    using TPtr      = typename T::Ptr;
    using Key       = std::string;

    // make sure T is a base class of Collection
    static_assert(std::is_base_of<Collection<ValueType>, T>::value, "T must be a base class of Collection");

    // all the keys
    Key  _all_key;

    // all the data
    TPtr _all;
    TPtr _cur;
    std::unordered_map<Key, TPtr> _data;

    bool has(const Key &name) {return has_key(name);}
    bool has_key(const Key &name) {return _data.find(name) != _data.end();}
    bool has_all() {return _all != nullptr;}
    bool has_any() {return _cur != nullptr;}
    TPtr all() {return _all;}
    TPtr get() {return _cur;}
    TPtr get(const Key &name) {return _data[name];}

    Sets(const Key &all_key="") : _all_key(all_key) {
        if (!_all_key.empty()) {
            _all = create(_all_key);
            _all->sorted(true);
            _all->duplicates(false);
        }
    }

    template <typename... Args>
    TPtr activate(const Key &name, Args... c) {
        Key key = name;
        if (key == "") {
            key = _all_key;
        }
        if (key == _all_key) {
            _cur = _all;
        } else {
            if (!has_key(key)) {
                _cur = create(name, c...);
            } else {
                _cur = _data[name];
            }
        }
        return _cur;
    }

    void add(const ValueType& item) {
        if (_cur) {
            _cur->add(item);
        }
        if (_all) {
            _all->add(item);
        }
    }

    // Add a range of values when ValueType is integral
    template <typename U = ValueType>
    typename std::enable_if<std::is_integral<U>::value>::type add(U first, U last, U step) {
        for (U value = first; value <= last; value += step) {
            add(value);
        }
    }

private:
    // create a new collection with additional values passed
    template <typename... Args>
    TPtr create(const Key &name, Args... c) {
        auto collection = std::make_shared<T>(name, c...);
        if (_all)
            collection->set_parent(_all);
        _data.emplace(name, collection);
        return collection;
    }
};

}    // namespace fem::model

#endif    // FEMASTER_SETS_H
