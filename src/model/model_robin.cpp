/**
 * @file model_robin.cpp
 * @brief Implements registration of Robin boundary conditions.
 */

#include "model.h"

#include "../core/logging.h"

namespace fem::model {

void Model::add_load(bc::Robin::Ptr load) {
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");
    logging::error(load != nullptr,
        "Model: cannot add a null Robin load");
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    _data->load_cols.get()->add(std::move(load));
}

} // namespace fem::model
