#pragma once
#include <string>
#include <unordered_map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "key_rule.h"

namespace fem::reader2 {

/**
 * \brief Fluent helper that validates keyword arguments against configured rules.
 */
struct KeyRules {
    // name -> rule (names are already uppercased by line parser)
    std::unordered_map<std::string, KeyRule> _rules;

    /// \brief Expose the underlying rule map for inspection.
    const std::unordered_map<std::string, KeyRule>& rules() const
    {
        return _rules;
    }

    /// \brief Provide descriptive text for documentation.
    KeyRules& describe(std::string name, std::string text)
    {
        _rules[std::move(name)].description = std::move(text);
        return *this;
    }

    /// \brief Mark a key as required.
    KeyRules& require(std::string name)
    {
        _rules[std::move(name)].required = true;
        return *this;
    }

    /// \brief Mark a key as optional and optionally supply a default.
    KeyRules& optional(std::string name, std::optional<std::string> def = std::nullopt)
    {
        auto& rule = _rules[std::move(name)];
        rule.required = false;
        rule.default_value = std::move(def);
        return *this;
    }

    /// \brief Restrict a key to a fixed set of uppercase values.
    KeyRules& allowed(std::string name, std::initializer_list<std::string> values)
    {
        auto& rule = _rules[name];
        rule.allowed_values.insert(values.begin(), values.end());
        return *this;
    }

    // Validate & materialize a final key map (string->string), applying defaults & enums.
    // Throws with a helpful message on violation.
    /// \brief Validate an input map and apply defaults.
    std::unordered_map<std::string, std::string>
    apply(const std::unordered_map<std::string,std::string>& input) const {
        std::unordered_map<std::string, std::string> out = input; // start with given keys
        for (const auto& [k, r] : _rules) {
            auto it = out.find(k);
            if (it == out.end()) {
                if (r.required) {
                    throw std::runtime_error("Missing required key: " + k);
                }
                if (r.default_value) out[k] = *r.default_value;
                continue;
            }
            if (!r.allowed_values.empty() && !r.allowed_values.count(it->second)) {
                std::ostringstream oss;
                oss << "Invalid value '" << it->second << "' for key " << k << ". Allowed: {";
                bool first=true; for (auto& v : r.allowed_values){ if(!first) oss<<","; first=false; oss<<v; } oss<<"}";
                throw std::runtime_error(oss.str());
            }
        }
        return out;
    }

    // Typed getter helper
    /// \brief Convert a string value into the requested type.
    template<class T>
    static T convert(const std::string& s){
        if constexpr (std::is_same_v<T,std::string>) return s;
        std::istringstream iss(s); T v{}; iss >> v;
        if (iss.fail()) throw std::runtime_error("Key conversion failed for value '"+s+"'");
        return v;
    }
};

} // namespace fem::reader2
