/******************************************************************************
 * @file data_storage.h
 * @brief Declares a sparse-like container storing FEM data by index.
 *
 * `DataStorage` behaves like an indexed array that lazily creates entries on
 * demand. It is used for nodal, elemental, and integration-point data where the
 * number of columns differs per field but the number of rows is fixed by the
 * model.
 *
 * @see src/data/data_storage.cpp
 * @see src/data/node_data_dict.h
 * @see src/data/elem_data_dict.h
 ******************************************************************************/

#pragma once

#include "../core/types_num.h"

#include <memory>
#include <vector>

namespace fem {
namespace model {

/******************************************************************************
 * @struct DataStorage
 * @brief Stores FEM matrices keyed by their row index.
 *
 * @tparam T Matrix-like type that provides `rows`, `cols`, `size`, `Zero`, and
 *           `(rows, cols)` constructors.
 ******************************************************************************/
template<typename T>
struct DataStorage {
    using Ptr = std::shared_ptr<DataStorage<T>>; ///< Shared pointer alias for storage blocks.
    using Key = Index;                           ///< Index type used for addressing entries.

    const Index _rows;        ///< Common row count shared by all stored matrices.
    std::vector<T> _data;     ///< Backing container storing the matrices.

    /******************************************************************************
     * @brief Constructs an empty storage with a prescribed row count.
     *
     * All slots are initialised with zero-sized matrices in order to delay
     * allocation until the first access.
     *
     * @param rows Number of rows every stored matrix must provide.
     ******************************************************************************/
    explicit DataStorage(Index rows)
        : _rows(rows) {
        _data.resize(rows, T::Zero(0, 0));
    }

    /******************************************************************************
     * @brief Erases the entry at `key` by resetting it to a zero-sized matrix.
     ******************************************************************************/
    void remove(Key key) {
        if (key >= _data.size()) {
            return;
        }
        _data[key] = T(0, 0);
    }

    /******************************************************************************
     * @brief Checks whether data has been created for `key`.
     ******************************************************************************/
    [[nodiscard]] bool has(Key key) const {
        if (key >= _data.size()) {
            return false;
        }
        if (_data[key].size() == 0) {
            return false;
        }
        return true;
    }

    /******************************************************************************
     * @brief Returns a mutable reference, creating the entry on demand.
     ******************************************************************************/
    T& get(Key key) {
        if (!has(key)) {
            reserve(key);
            _data[key] = T(0, 0);
        }
        return _data[key];
    }

    /******************************************************************************
     * @brief Returns a copy of the stored data or an empty matrix when missing.
     ******************************************************************************/
    T get(Key key) const {
        if (!has(key)) {
            return T(0, 0);
        }
        return _data[key];
    }

    /// Returns a mutable reference to the entry via array syntax.
    T& operator[](Key key) { return get(key); }

    /// Returns a copy of the entry via array syntax (const overload).
    T operator[](Key key) const { return get(key); }

    /******************************************************************************
     * @brief Allocates a matrix with `_rows` rows and `entries` columns.
     *
     * Existing data is left untouched.
     ******************************************************************************/
    void create(Key key, int entries) {
        if (!has(key)) {
            reserve(key);
            _data[key] = T(_rows, entries);
        }
    }

private:
    /******************************************************************************
     * @brief Resizes the storage to accommodate the provided key.
     ******************************************************************************/
    void reserve(Key key) {
        if (key >= _data.size()) {
            _data.resize(key + 1, T::Zero(0, 0));
        }
    }
};

} // namespace model
} // namespace fem

