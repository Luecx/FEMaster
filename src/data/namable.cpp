/******************************************************************************
 * @file namable.cpp
 * @brief Implements the `Namable` mixin.
 *
 * Contains the trivial constructor definition required by multiple translation
 * units to avoid duplicate symbol issues when linking.
 *
 * @see src/data/namable.h
 * @see src/data/collection.h
 ******************************************************************************/

#include "namable.h"

#include <utility>

namespace fem {
namespace model {

/******************************************************************************
 * @copydoc Namable::Namable
 ******************************************************************************/
Namable::Namable(std::string p_name)
    : name(std::move(p_name)) {}

} // namespace model
} // namespace fem

