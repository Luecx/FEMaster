#pragma once

#include <optional>
#include <string>
#include <unordered_set>

namespace fem::reader2 {

/**
 * \brief Describes validation requirements for a single keyword argument.
 */
struct KeyRule {
    bool                              required = false;           ///< Whether the key must appear.
    std::optional<std::string>        default_value;              ///< Default value if optional.
    std::unordered_set<std::string>   allowed_values;             ///< Enumerated allowed values.
    std::string                       description;                ///< Human-facing documentation.
};

} // namespace fem::reader2
