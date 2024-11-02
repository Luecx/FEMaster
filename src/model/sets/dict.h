//
// Created by Finn Eggers on 02.11.24.
//

#ifndef FEMASTER_DICT_H
#define FEMASTER_DICT_H

#include <type_traits>
#include <memory>
#include <unordered_map>

namespace fem::model {

template <typename T>
struct Dict {
    using Key = std::string;
    using Ptr = std::shared_ptr<T>;

    // Data map holding instances of T by key
    std::unordered_map<Key, Ptr> _data;
    Ptr _cur;

    bool has(const Key &name) { return has_key(name); }
    bool has_key(const Key &name) { return _data.find(name) != _data.end(); }
    bool has_any() { return _cur != nullptr; }
    Ptr get() { return _cur; }
    Ptr get(const Key &name) { return _data.at(name); }

    template <typename Derived = T, typename... Args>
    Ptr activate(const Key &name, Args... args) {
        static_assert(std::is_base_of<T, Derived>::value, "Derived must be derived from T");

        if (!has_key(name)) {
            _cur = create<Derived>(name, args...);
        } else {
            _cur = _data[name];
        }
        return _cur;
    }

    private:
    // Create a new instance of the specified type with additional values passed
    template <typename Derived = T, typename... Args>
    Ptr create(const Key &name, Args... args) {
        static_assert(std::is_base_of<T, Derived>::value, "Derived must be derived from T");

        auto instance = std::make_shared<Derived>(args...);
        _data.emplace(name, instance);
        return instance;
    }
};

}    // namespace fem::model

#endif    // FEMASTER_DICT_H
