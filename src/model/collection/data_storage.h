//
// Created by f_eggers on 08.11.2024.
//

#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H

#include "../../core/core.h"

namespace fem::model {

template<typename T>
struct DataStorage {
    using Ptr = std::shared_ptr<DataStorage<T>>;
    using Key = Index;

    std::vector<T> _data;
       
    void remove(Index key) {
        if (key >= _data.size()) {
            return;
        }
        _data[key] = T(0, 0);
    }
       
    bool has(Index key) const {
        if (key >= _data.size()) {
            return false;
        }
        if (_data[key].size() == 0) {
            return false;
        }
        return true;
    }

    T& get(Index key) {
        if (!has(key)) {
            reserve(key);
            _data[key] = T(0, 0);
        }
        return _data[key];
    }
    T get(Index key) const {
        if (!has(key)) {
            return T(0, 0);
        }
        return _data[key];
    }

    T& operator[](Index key) {
        return get(key);
    }
    T operator[](Index key) const {
        return get(key);
    }

    private:
    void reserve(Index key) {
        if (key >= _data.size()) {
            _data.resize(key + 1);
        }
    }
};

}


#endif //DATA_STORAGE_H
