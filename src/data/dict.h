#ifndef FEMASTER_DICT_H
#define FEMASTER_DICT_H

#include <type_traits>
#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include "namable.h"

namespace fem::model {

template <typename T, typename Key = std::string>
struct Dict {
    using Ptr  = std::shared_ptr<Dict<T, Key>>;
    using TPtr = std::shared_ptr<T>;

    // Data map/list holding instances of T by key
    std::conditional_t<std::is_same<Key, std::string>::value,
                       std::unordered_map<Key, TPtr>,
                       std::vector<TPtr>> _data;

    // if Key is string, make sure the entity is namable
    static_assert(!std::is_same_v<Key, std::string> || std::is_base_of_v<Namable, T>,
                      "Key must be a string if T is namable");
    TPtr _cur;

    bool has(const Key &key) {
        if constexpr (std::is_same<Key, std::string>::value) {
            return has_key(key);
        } else {
            return key < _data.size() && _data[key] != nullptr;
        }
    }

    bool has_key(const Key &key) {
        if constexpr (std::is_same<Key, std::string>::value) {
            return _data.find(key) != _data.end();
        }
        return false;
    }

    bool has_any() { return _cur != nullptr; }

    TPtr get() { return _cur; }

    TPtr get(const Key &key) {
        if (!has(key))
            return nullptr;

        if constexpr (std::is_same<Key, std::string>::value) {
            return _data.at(key);
        } else {
            return key < _data.size() ? _data[key] : nullptr;
        }
    }

    template <typename Derived = T, typename... Args>
    TPtr activate(const Key &key, Args... args) {
        static_assert(std::is_base_of<T, Derived>::value, "Derived must be derived from T");

        if (!has(key)) {
            _cur = create<Derived>(key, args...);
        } else {
            _cur = get(key);
        }
        return _cur;
    }

    void remove(const Key &key) {
        if constexpr (std::is_same<Key, std::string>::value) {
            _data.erase(key);
        } else {
            if (key < _data.size()) _data[key] = nullptr;
        }
    }

    // Create a new instance and store it by key
    template <typename Derived = T, typename... Args>
    TPtr create(const Key &key, Args... args) {
        static_assert(std::is_base_of<T, Derived>::value, "Derived must be derived from T");

        if constexpr (std::is_same<Key, std::string>::value) {
            auto instance = std::make_shared<Derived>(key, args...);
            _data.emplace(key, instance);
            return instance;
        } else {
            auto instance = std::make_shared<Derived>(args...);
            if (key >= _data.size()) _data.resize(key + 1);
            _data[key] = instance;
            return instance;
        }
    }

    // begin and end function
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }
};

} // namespace fem::model

#endif // FEMASTER_DICT_H
