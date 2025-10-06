#pragma once
/**
 * @file key_rules.h
 * @brief Declares KeyRules, a fluent helper to validate keyword arguments
 *        against configured rules (required/optional, defaults, enums).
 *
 * @date 06.10.2025
 * @version 1.0
 */

#include <string>
#include <unordered_map>
#include <optional>
#include <initializer_list>
#include <sstream>
#include <type_traits>

#include "key_rule.h"

namespace fem::reader2 {

/**
 * @brief Fluent helper that validates keyword arguments against configured rules.
 *
 * Names are assumed to be uppercased by the line parser before validation.
 */
struct KeyRules {
    /// @brief Underlying rule map (name → rule).
    std::unordered_map<std::string, KeyRule> _rules;

    /** @brief Expose the underlying rule map for inspection. */
    const std::unordered_map<std::string, KeyRule>& rules() const;

    /** @brief Provide descriptive text for documentation. */
    KeyRules& describe(std::string name, std::string text);

    /** @brief Mark a key as required. */
    KeyRules& require(std::string name);

    /**
     * @brief Mark a key as optional and optionally supply a default.
     * @param name Key name
     * @param def  Optional default value (if not provided by the user)
     */
    KeyRules& optional(std::string name, std::optional<std::string> def = std::nullopt);

    /**
     * @brief Restrict a key to a fixed set of uppercase values.
     * @param name   Key name
     * @param values Allowed uppercase values
     */
    KeyRules& allowed(std::string name, std::initializer_list<std::string> values);

    /**
     * @brief Validate an input key→value map and apply defaults/enums.
     * @throws std::runtime_error on violation (missing required, invalid enum, etc.).
     */
    std::unordered_map<std::string, std::string>
    apply(const std::unordered_map<std::string, std::string>& input) const;

    /**
     * @brief Convert a string value into the requested type.
     * @tparam T Target type; `std::string` is returned verbatim, arithmetic types via stream.
     * @throws std::runtime_error if conversion fails.
     */
    template<class T>
    static T convert(const std::string& s) {
        if constexpr (std::is_same_v<T, std::string>) {
            return s;
        } else {
            std::istringstream iss(s);
            T v{};
            iss >> v;
            if (iss.fail()) {
                throw std::runtime_error("Key conversion failed for value '" + s + "'");
            }
            return v;
        }
    }

    /** @brief Factory helper. */
    static KeyRules make();
};

} // namespace fem::reader2
