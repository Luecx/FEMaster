/**
 * @file key_rules.cpp
 * @brief Implements KeyRules for keyword validation and defaults application.
 *
 * @date 06.10.2025
 * @version 1.0
 */

#include "key_rules.h"
#include <stdexcept>
#include <sstream>

namespace fem::reader2 {

const std::unordered_map<std::string, KeyRule>& KeyRules::rules() const {
    return _rules;
}

KeyRules& KeyRules::describe(std::string name, std::string text) {
    _rules[std::move(name)].description = std::move(text);
    return *this;
}

KeyRules& KeyRules::require(std::string name) {
    _rules[std::move(name)].required = true;
    return *this;
}

KeyRules& KeyRules::optional(std::string name, std::optional<std::string> def) {
    auto& rule = _rules[std::move(name)];
    rule.required = false;
    rule.default_value = std::move(def);
    return *this;
}

KeyRules& KeyRules::allowed(std::string name, std::initializer_list<std::string> values) {
    auto& rule = _rules[name];
    rule.allowed_values.insert(values.begin(), values.end());
    return *this;
}

std::unordered_map<std::string, std::string>
KeyRules::apply(const std::unordered_map<std::string, std::string>& input) const {
    std::unordered_map<std::string, std::string> out = input; // start with given keys
    for (const auto& [k, r] : _rules) {
        auto it = out.find(k);
        if (it == out.end()) {
            if (r.required) {
                throw std::runtime_error("Missing required key: " + k);
            }
            if (r.default_value) {
                out[k] = *r.default_value;
            }
            continue;
        }
        if (!r.allowed_values.empty() && !r.allowed_values.count(it->second)) {
            std::ostringstream oss;
            oss << "Invalid value '" << it->second << "' for key " << k << ". Allowed: {";
            bool first = true;
            for (const auto& v : r.allowed_values) {
                if (!first) oss << ",";
                first = false;
                oss << v;
            }
            oss << "}";
            throw std::runtime_error(oss.str());
        }
    }
    return out;
}

KeyRules KeyRules::make() {
    return {};
}

} // namespace fem::reader2
